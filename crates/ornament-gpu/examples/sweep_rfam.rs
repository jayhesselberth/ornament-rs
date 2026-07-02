//! All-Rfam MSV sweep: validate the batch MSV kernel across the *real* Rfam model-length (M)
//! distribution, not just tRNA. For every model in an Infernal `.cm`/`.cm.gz` collection we build
//! the uint8 MSV profile from its embedded fp7 filter, score a batch of windows over a resident
//! genome strand on the GPU, and record throughput + which kernel ran (shared-mem vs global
//! fallback). A sample of models is checked for parity against the CPU striped filter.
//!
//! This answers: does the shared-memory speedup hold across Rfam, and how much of the collection
//! falls past the 48 KB shared-mem crossover into the (slow) global-fallback / warp-per-window
//! regime? Run on a `compgpu` node:
//!
//! ```text
//! RFAM_CM=~/.cache/rfam_models/Rfam.cm.gz \
//!   cargo run --release --example sweep_rfam --features cuda
//! ```
//!
//! Env: `RFAM_CM` (path to the .cm/.cm.gz, required), `BENCH_STRAND` (genome residues, default
//! 8M — smaller than bench_msv since we scan it once per model), `BENCH_W` (window, default 100),
//! `SWEEP_PARITY` (# models to parity-check vs CPU, default 40), `SWEEP_MAX_MODELS` (cap, default
//! all). `ORNAMENT_GPU_SMEM=0` forces the global kernel for an A/B of the whole collection.

fn env_usize(name: &str, def: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(def)
}

#[cfg(feature = "cuda")]
fn main() {
    use ornament_alphabet::Alphabet;
    use ornament_gpu::{
        nats_to_bits, ByteMsvProfile, DeviceStrand, FlatProfile, GpuContext, Tiles,
    };
    use ornament_hmm::profile::{bg_freqs, P7Profile};
    use ornament_hmm::MsvProfile;
    use std::time::Instant;

    let cm_path = std::env::var("RFAM_CM").unwrap_or_else(|_| {
        eprintln!("set RFAM_CM to an Infernal .cm/.cm.gz (e.g. ~/.cache/rfam_models/Rfam.cm.gz)");
        std::process::exit(1);
    });
    let strand_len = env_usize("BENCH_STRAND", 8_000_000);
    let w = env_usize("BENCH_W", 100);
    let parity_n = env_usize("SWEEP_PARITY", 40);
    let max_models = env_usize("SWEEP_MAX_MODELS", usize::MAX);
    let smem_off = std::env::var("ORNAMENT_GPU_SMEM").as_deref() == Ok("0");
    // Mirror the kernel's real dispatch: the context opts into the sm_80 ~164 KB per-block shared
    // limit and picks the largest block size {128,64,32} whose (Kp+bd)*(M+1) fits — so a model runs
    // on shared memory unless even 32 threads overflow ~164 KB. (166912 = A30 opt-in max.)
    let smem_fits = |m: usize, kp: usize| {
        [128usize, 64, 32]
            .iter()
            .any(|&bd| (kp + bd) * (m + 1) <= 166_912)
    };

    match ornament_gpu::device_count() {
        Ok(0) | Err(_) => {
            eprintln!("no CUDA device visible — run on a compgpu node");
            std::process::exit(1);
        }
        Ok(c) => eprintln!("CUDA devices: {c}"),
    }

    // Parse every model in the collection (handles gzip + multi-record).
    eprintln!("parsing {cm_path} ...");
    let models = ornament_scfg::parse_cm_records_file(&cm_path).unwrap_or_else(|e| {
        eprintln!("parse error: {e}");
        std::process::exit(1);
    });
    eprintln!("{} models", models.len());

    let abc = Alphabet::rna();

    // One resident genome strand, scored by every model.
    let mut s: u64 = 0x243F6A8885A308D3;
    let bases = [b'A', b'C', b'G', b'U'];
    let seq: String = (0..strand_len)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            bases[((s >> 40) & 3) as usize] as char
        })
        .collect();
    let strand = abc.digitize(&seq).expect("digitize");
    let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
    // One reusable context for every model — hoists the stream/pinned-buffer setup out of the
    // per-model path (that setup dominated the earlier per-call sweep).
    let ctx = GpuContext::new().expect("gpu context");

    // Tile the strand once (length 2W, step W) — reused for every model.
    let tile = (2 * w).max(1);
    let step = w.max(1);
    let mut spans: Vec<(usize, usize)> = Vec::new();
    let mut start = 1usize;
    loop {
        let end = (start + tile - 1).min(strand_len);
        spans.push((start, end));
        if end >= strand_len {
            break;
        }
        start += step;
    }
    let tiles = Tiles::from_spans(&spans);
    let n_tiles = tiles.len();
    eprintln!("strand={strand_len}  W={w}  tiles/model={n_tiles}  smem_forced_off={smem_off}");

    // Buckets over the M distribution (aligns with the shared-mem crossover at M≈336).
    let bucket = |m: usize| match m {
        0..=99 => 0,
        100..=199 => 1,
        200..=335 => 2,
        336..=999 => 3,
        _ => 4,
    };
    let bnames = ["M<100", "M100-199", "M200-335", "M336-999", "M>=1000"];
    let mut b_count = [0usize; 5];
    let mut b_secs = [0.0f64; 5];
    let mut b_tiles = [0.0f64; 5];

    let mut ok_models = 0usize;
    let mut skipped = 0usize;
    let mut fits_smem = 0usize;
    let mut falls_global = 0usize;
    let mut parity_checked = 0usize;
    let mut parity_fail = 0usize;
    let mut max_m = 0usize;
    let total_t0 = Instant::now();

    for (mi, cm) in models.iter().enumerate().take(max_models) {
        let Some(fp7) = cm.fp7_text.as_deref() else {
            skipped += 1;
            continue;
        };
        let hmm = match ornament_hmm::parse_p7_hmm(fp7) {
            Ok(h) => h,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };
        let bg = bg_freqs(hmm.k);
        let prof = P7Profile::config_local(&hmm, &abc, &bg, 2 * w);
        let bp = ByteMsvProfile::new(&prof);
        let m = bp.m;
        max_m = max_m.max(m);
        if smem_off || !smem_fits(m, bp.kp) {
            falls_global += 1;
        } else {
            fits_smem += 1;
        }

        // One timed batch per model (the per-call stream/pinned setup makes extra reps costly at
        // 4k models); the first model's call also warms up the kernel JIT.
        let t0 = Instant::now();
        let gpu = match ctx.msv_u8_nats(&bp, &dstrand, &tiles) {
            Ok(o) => o,
            Err(e) => {
                eprintln!("model {} ({}): gpu error {e}", mi, cm.name);
                skipped += 1;
                continue;
            }
        };
        let best = t0.elapsed().as_secs_f64();
        ok_models += 1;
        let b = bucket(m);
        b_count[b] += 1;
        b_secs[b] += best;
        b_tiles[b] += n_tiles as f64;

        // Parity-check a spread of models against the CPU striped filter (reuse the batch above;
        // first tiles only, cheap).
        if parity_checked < parity_n && mi % (models.len() / parity_n.max(1)).max(1) == 0 {
            let mut sp = MsvProfile::new(&prof);
            sp.set_length(2 * w);
            let _ = FlatProfile::new(&prof); // (exercise the f32 builder path too)
            let sample = 200.min(spans.len());
            let mut fails = 0;
            for k in 0..sample {
                let (a, e) = spans[k];
                let mut win = Vec::with_capacity(e - a + 3);
                win.push(ornament_alphabet::alphabet::SENTINEL);
                win.extend_from_slice(&strand[a..=e]);
                win.push(ornament_alphabet::alphabet::SENTINEL);
                let cpu_bits = sp.score_bits(&win);
                let gpu_bits = nats_to_bits(gpu[k], e - a + 1);
                let agree = (gpu_bits - cpu_bits).abs() < 0.01
                    || (gpu_bits.is_infinite() && cpu_bits.is_infinite());
                if !agree {
                    fails += 1;
                }
            }
            parity_checked += 1;
            if fails > 0 {
                parity_fail += 1;
                eprintln!(
                    "PARITY FAIL model {} ({}) M={m}: {fails}/{sample}",
                    mi, cm.name
                );
            }
        }
    }

    let total_secs = total_t0.elapsed().as_secs_f64();
    println!("\n=== all-Rfam MSV sweep ===");
    println!(
        "models: {ok_models} scored, {skipped} skipped, max M = {max_m}  (smem_forced_off={smem_off})"
    );
    println!("kernel dispatch: {fits_smem} fit shared-mem, {falls_global} fell to global");
    println!("parity: {parity_checked} models checked, {parity_fail} failed");
    println!("\nthroughput by model length (M):");
    println!(
        "  {:<10} {:>7} {:>12} {:>16}",
        "bucket", "models", "sum secs", "M tiles/s"
    );
    for b in 0..5 {
        if b_count[b] == 0 {
            continue;
        }
        let tps = if b_secs[b] > 0.0 {
            b_tiles[b] / b_secs[b] / 1e6
        } else {
            0.0
        };
        println!(
            "  {:<10} {:>7} {:>12.4} {:>16.2}",
            bnames[b], b_count[b], b_secs[b], tps
        );
    }
    let all_tiles: f64 = b_tiles.iter().sum();
    let all_secs: f64 = b_secs.iter().sum();
    println!(
        "\ntotal MSV kernel time {:.2}s (+parse/setup = {:.1}s wall); overall {:.1} M tiles/s",
        all_secs,
        total_secs,
        all_tiles / all_secs / 1e6
    );
}

#[cfg(not(feature = "cuda"))]
fn main() {
    let _ = env_usize;
    eprintln!("sweep_rfam needs --features cuda (GPU path compiled out)");
}
