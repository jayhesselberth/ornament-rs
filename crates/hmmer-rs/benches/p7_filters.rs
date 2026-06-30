//! Microbenchmarks for the p7 HMM filter DP kernels (MSV / Viterbi / Forward).
//!
//! These are the hot paths the SIMD/rayon work (#5) will accelerate; this establishes the
//! scalar baseline and a regression guard. The model is the RF00005 filter HMM (M=71); the
//! target is a deterministic ~3 kb RNA sequence so timings reflect realistic throughput.
//!
//! Run:
//!   srun -p rna -c 24 cargo bench -p hmmer-rs --bench p7_filters
//!   cargo bench -p hmmer-rs --bench p7_filters -- forward

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use easel_rs::Alphabet;
use hmmer_rs::profile::bg_freqs;
use hmmer_rs::{
    forward_bits, forward_nats, msv_bits, parse_p7_hmm, viterbi_bits, P7Hmm, P7Profile,
};

/// Deterministic pseudo-random RNA text of length `n` (xorshift; no external RNG).
fn pseudo_rna(n: usize, seed: u64) -> String {
    let mut s = seed | 1;
    let bases = [b'A', b'C', b'G', b'U'];
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            bases[(s % 4) as usize] as char
        })
        .collect()
}

fn load_hmm() -> P7Hmm {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/benches/data/trna_fp7.hmm"
    ))
    .expect("read filter HMM");
    parse_p7_hmm(&text).expect("parse filter HMM")
}

fn bench_filters(c: &mut Criterion) {
    let hmm = load_hmm();
    let abc = Alphabet::rna();
    let l = 3000;
    let dsq = abc.digitize(&pseudo_rna(l, 0x5eed)).expect("digitize");

    let mut g = c.benchmark_group("p7_filters");
    g.throughput(Throughput::Elements(l as u64)); // residues/sec

    g.bench_function("msv", |b| {
        b.iter(|| msv_bits(black_box(&hmm), black_box(&abc), black_box(&dsq)))
    });
    g.bench_function("viterbi", |b| {
        b.iter(|| viterbi_bits(black_box(&hmm), black_box(&abc), black_box(&dsq)))
    });
    g.bench_function("forward", |b| {
        b.iter(|| forward_bits(black_box(&hmm), black_box(&abc), black_box(&dsq)))
    });
    g.finish();

    // DP-only comparison (profile pre-built, as in the windowing pipeline): the log-space
    // Forward vs the striped-SSE odds-space Forward, isolating the SIMD kernel speedup.
    let prof = P7Profile::config_local(&hmm, &abc, &bg_freqs(hmm.k), l);
    let mut gd = c.benchmark_group("p7_forward_dp");
    gd.throughput(Throughput::Elements(l as u64));
    gd.bench_function("log_space", |b| {
        b.iter(|| forward_nats(black_box(&prof), black_box(&dsq)))
    });
    #[cfg(target_arch = "x86_64")]
    {
        let sp = hmmer_rs::StripedProfile::new(&prof);
        gd.bench_function("striped_sse", |b| b.iter(|| sp.score_nats(black_box(&dsq))));
    }
    gd.finish();
}

criterion_group!(benches, bench_filters);
criterion_main!(benches);
