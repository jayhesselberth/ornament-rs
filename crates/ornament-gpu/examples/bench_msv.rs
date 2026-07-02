//! Throughput benchmark for the batch MSV filter, used to measure the streaming/double-buffer
//! speedup and the GPU-vs-CPU ratio. Build/run on a `compgpu` node:
//!
//! ```text
//! # A/B the stream overlap (1 = no overlap baseline, 2 = double-buffered pipeline):
//! ORNAMENT_GPU_STREAMS=1 cargo run --release --example bench_msv --features cuda
//! ORNAMENT_GPU_STREAMS=2 cargo run --release --example bench_msv --features cuda
//! ```
//!
//! Env knobs: `ORNAMENT_GPU_STREAMS` (1..8), `ORNAMENT_GPU_CHUNK` (tiles/chunk),
//! `BENCH_STRAND` (strand residues, default 64M), `BENCH_W` (window/step, default 100),
//! `BENCH_REPS` (timed reps, default 8).
//!
//! Without `--features cuda` it prints a notice and exits (the GPU path is compiled out).

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
        nats_to_bits, ByteMsvProfile, DeviceStrand, FlatProfile, GpuContext, GpuError, Tiles,
        VitFlatProfile,
    };
    use ornament_hmm::profile::{bg_freqs, P7Profile};
    use std::time::Instant;

    let strand_len = env_usize("BENCH_STRAND", 64_000_000);
    let w = env_usize("BENCH_W", 100);
    let reps = env_usize("BENCH_REPS", 8);
    // Which kernel to bench: "u8" (default, reduced-precision MSV), "f32" (MSV), or "vit" (Viterbi).
    let kernel = std::env::var("BENCH_KERNEL").unwrap_or_else(|_| "u8".into());
    let streams = std::env::var("ORNAMENT_GPU_STREAMS").unwrap_or_else(|_| "2 (default)".into());
    let chunk = std::env::var("ORNAMENT_GPU_CHUNK").unwrap_or_else(|_| "32768 (default)".into());

    // Model: the embedded tRNA filter HMM (M≈71) — the small-model case the scan stresses most.
    let hmm_path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../ornament-scfg/tests/data/trna_fp7.hmm"
    );
    let hmm = ornament_hmm::parse_p7_hmm(&std::fs::read_to_string(hmm_path).expect("read hmm"))
        .expect("parse hmm");
    let abc = Alphabet::rna();
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(&hmm, &abc, &bg, 100);
    let bp = ByteMsvProfile::new(&prof);
    let fp = FlatProfile::new(&prof);
    let vp = VitFlatProfile::new(&prof);

    // Random strand.
    let mut s: u64 = 0x9E3779B97F4A7C15;
    let bases = [b'A', b'C', b'G', b'U'];
    let seq: String = (0..strand_len)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            bases[((s >> 40) & 3) as usize] as char
        })
        .collect();
    let strand = abc.digitize(&seq).expect("digitize");

    // Tile like the pipeline: length 2W, step W (2× overlap), 1-based inclusive spans.
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
    let n = tiles.len();

    match ornament_gpu::device_count() {
        Ok(0) | Err(_) => {
            eprintln!("no CUDA device visible — run on a compgpu node");
            std::process::exit(1);
        }
        Ok(c) => eprintln!("CUDA devices: {c}"),
    }

    println!(
        "strand={strand_len} residues  W={w}  tiles={n}  streams={streams}  chunk={chunk}  reps={reps}"
    );
    // Counter-free occupancy/register report (ncu Occupancy substitute where perf counters are
    // locked). BENCH_INFO=1 prints it for this model's (Kp, M).
    if std::env::var("BENCH_INFO").as_deref() == Ok("1") {
        ornament_gpu::print_kernel_info(prof.kp, prof.m);
    }

    let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
    // Reusable context (created once) — realistic many-batch perf and a stable target for nsys/ncu.
    let ctx = GpuContext::new().expect("gpu context");

    // One scoring call for the selected kernel.
    let score = |ctx: &GpuContext| -> Result<Vec<f32>, GpuError> {
        match kernel.as_str() {
            "f32" => ctx.msv_nats(&fp, &dstrand, &tiles),
            "vit" => ctx.vit_nats(&vp, &dstrand, &tiles),
            _ => ctx.msv_u8_nats(&bp, &dstrand, &tiles),
        }
    };

    // Warm up (kernel JIT / first-touch allocations) and sanity-check the result is finite-or-inf.
    let warm = score(&ctx).expect("gpu warmup");
    let _ = nats_to_bits(warm[0], spans[0].1 - spans[0].0 + 1);

    let mut best = f64::INFINITY;
    let mut sum = 0.0f64;
    for _ in 0..reps {
        let t0 = Instant::now();
        let out = score(&ctx).unwrap_or_else(|e: GpuError| {
            eprintln!("gpu error: {e}");
            std::process::exit(1);
        });
        std::hint::black_box(&out);
        let dt = t0.elapsed().as_secs_f64();
        best = best.min(dt);
        sum += dt;
    }
    let avg = sum / reps as f64;
    let tiles_per_s = n as f64 / best;
    let cells = n as f64 * (2 * w) as f64 * prof.m as f64; // DP cells scored (tiles × L × M)

    println!("GPU kernel={kernel}:");
    println!("  best {:.4}s   avg {:.4}s", best, avg);
    println!(
        "  {:.2} M tiles/s   {:.2} G DP-cells/s",
        tiles_per_s / 1e6,
        cells / best / 1e9
    );

    // Single-core CPU reference on a sample (the striped uint8 filter), extrapolated. The real
    // pipeline runs this under rayon, so multiply by core count for the practical CPU number.
    #[cfg(target_arch = "x86_64")]
    {
        use ornament_hmm::MsvProfile;
        let sample = n.min(20_000);
        let mut sp = MsvProfile::new(&prof);
        sp.set_length(2 * w);
        let t0 = Instant::now();
        let mut acc = 0.0f32;
        for &(a, b) in spans.iter().take(sample) {
            let mut win = Vec::with_capacity(b - a + 3);
            win.push(ornament_alphabet::alphabet::SENTINEL);
            win.extend_from_slice(&strand[a..=b]);
            win.push(ornament_alphabet::alphabet::SENTINEL);
            acc += sp.score_bits(&win);
        }
        std::hint::black_box(acc);
        let cpu_dt = t0.elapsed().as_secs_f64();
        let cpu_tiles_per_s = sample as f64 / cpu_dt;
        println!("CPU u8 MSV (1 core, striped SIMD), {sample} tiles:");
        println!(
            "  {:.3} M tiles/s   GPU is {:.1}× one core",
            cpu_tiles_per_s / 1e6,
            tiles_per_s / cpu_tiles_per_s
        );
    }
}

#[cfg(not(feature = "cuda"))]
fn main() {
    let _ = env_usize;
    eprintln!("bench_msv needs --features cuda (GPU path compiled out)");
}
