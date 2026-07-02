//! Tile-aggregation benchmark: score many short sequences against one model, comparing
//! **K separate GPU launches** (one per sequence) vs **one aggregated launch** over a concatenated
//! resident strand (`DeviceStrand::upload_concat` + `Tiles::push_segment`). ncu showed a single
//! small batch underfills the A30 (few blocks, ~0.2 waves); aggregation packs many sequences'
//! tiles into one launch to fill it. Run on a `compgpu` node:
//!
//! ```text
//! cargo run --release --example bench_agg --features cuda
//! ```
//! Env: `AGG_SEGMENTS` (# sequences, default 256), `AGG_SEGLEN` (residues each, default 4000),
//! `BENCH_W` (window, default 100), `BENCH_REPS` (default 5).

fn env_usize(name: &str, def: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(def)
}

#[cfg(feature = "cuda")]
fn main() {
    use ornament_alphabet::Alphabet;
    use ornament_gpu::{ByteMsvProfile, DeviceStrand, GpuContext, Tiles};
    use ornament_hmm::profile::{bg_freqs, P7Profile};
    use std::time::Instant;

    let n_seg = env_usize("AGG_SEGMENTS", 256);
    let seg_len = env_usize("AGG_SEGLEN", 4000);
    let w = env_usize("BENCH_W", 100);
    let reps = env_usize("BENCH_REPS", 5);

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

    match ornament_gpu::device_count() {
        Ok(0) | Err(_) => {
            eprintln!("no CUDA device visible — run on a compgpu node");
            std::process::exit(1);
        }
        Ok(c) => eprintln!("CUDA devices: {c}"),
    }

    // n_seg short random sequences (sentinel-padded), each tiled length-2W step-W.
    let mut s: u64 = 0x2545F4914F6CDD1D;
    let bases = [b'A', b'C', b'G', b'U'];
    let seqs: Vec<Vec<u8>> = (0..n_seg)
        .map(|_| {
            let text: String = (0..seg_len)
                .map(|_| {
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                    bases[((s >> 40) & 3) as usize] as char
                })
                .collect();
            abc.digitize(&text).expect("digitize")
        })
        .collect();
    let spans: Vec<(usize, usize)> = {
        let tile = (2 * w).max(1);
        let step = w.max(1);
        let mut v = Vec::new();
        let mut start = 1usize;
        loop {
            let end = (start + tile - 1).min(seg_len);
            v.push((start, end));
            if end >= seg_len {
                break;
            }
            start += step;
        }
        v
    };
    let tiles_total = n_seg * spans.len();
    println!("segments={n_seg} seglen={seg_len} W={w} tiles/seg={} total_tiles={tiles_total} reps={reps}", spans.len());

    let ctx = GpuContext::new().expect("ctx");

    // (a) SEPARATE: one upload + one launch per sequence (the underfilled pattern).
    let sep = |ctx: &GpuContext| {
        for seq in &seqs {
            let ds = DeviceStrand::upload(seq).expect("upload");
            let t = Tiles::from_spans(&spans);
            std::hint::black_box(ctx.msv_u8_nats(&bp, &ds, &t).expect("sep"));
        }
    };
    // (b) AGGREGATED: one concatenated upload + one launch over all sequences' tiles.
    let refs: Vec<&[u8]> = seqs.iter().map(|v| v.as_slice()).collect();
    let agg = |ctx: &GpuContext| {
        let (ds, bases) = DeviceStrand::upload_concat(&refs).expect("concat");
        let mut t = Tiles::default();
        for &b in &bases {
            t.push_segment(b, &spans);
        }
        std::hint::black_box(ctx.msv_u8_nats(&bp, &ds, &t).expect("agg"));
    };

    sep(&ctx);
    agg(&ctx); // warm up both paths

    let time = |f: &dyn Fn(&GpuContext)| {
        let mut best = f64::INFINITY;
        for _ in 0..reps {
            let t0 = Instant::now();
            f(&ctx);
            best = best.min(t0.elapsed().as_secs_f64());
        }
        best
    };
    let t_sep = time(&sep);
    let t_agg = time(&agg);

    println!(
        "separate ({n_seg} launches): {:.4}s  {:.2} M tiles/s",
        t_sep,
        tiles_total as f64 / t_sep / 1e6
    );
    println!(
        "aggregated (1 launch):       {:.4}s  {:.2} M tiles/s",
        t_agg,
        tiles_total as f64 / t_agg / 1e6
    );
    println!("speedup: {:.1}x", t_sep / t_agg);
}

#[cfg(not(feature = "cuda"))]
fn main() {
    let _ = env_usize;
    eprintln!("bench_agg needs --features cuda (GPU path compiled out)");
}
