//! Microbenchmarks for the covariance-model DP kernels — the hot paths the SIMD/rayon
//! acceleration (#5) will target. Establishes the scalar baseline + a regression guard.
//!
//! - `cm_dp`: the scanning CYK / Inside DP and the CYK alignment traceback (the inner
//!   recursions that dominate runtime), on the embedded-tRNA fixture.
//! - `cm_search`: the realistic both-strand multi-hit search vs. the p7-filtered pipeline,
//!   on a longer synthetic sequence, to show the filter's pruning benefit.
//!
//! Run:
//!   srun -p rna -c 24 cargo bench -p infernal-rs --bench cm_search
//!   cargo bench -p infernal-rs --bench cm_search -- cyk

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use easel_rs::{read_fasta, Alphabet, Dsq};
use infernal_rs::{
    align_glocal, cm_pipeline_search, configure_local, configure_scores, cyk_scan, cyk_search,
    inside_scan, parse_cm_file, Cm, EmitMap, PipelineParams,
};

fn data(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn load_cm() -> Cm {
    parse_cm_file(data("tRNA.cm")).expect("parse tRNA.cm")
}

fn read_seq(abc: &Alphabet, file: &str) -> Vec<Dsq> {
    let recs = read_fasta(data(file)).expect("read fasta");
    abc.digitize(&recs[0].seq).expect("digitize")
}

/// Deterministic ~`n`-nt RNA with the tRNA consensus embedded at two positions.
fn synthetic_with_trnas(abc: &Alphabet, consensus: &str, n: usize) -> Vec<Dsq> {
    let mut s: u64 = 0xC0FFEE;
    let bases = [b'A', b'C', b'G', b'U'];
    let mut seq: String = (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            bases[(s % 4) as usize] as char
        })
        .collect();
    // Splice the consensus in at ~1/4 and ~3/4.
    let q = n / 4;
    seq.replace_range(q..q + consensus.len(), consensus);
    seq.replace_range(3 * q..3 * q + consensus.len(), consensus);
    abc.digitize(&seq).expect("digitize")
}

fn bench_cm_dp(c: &mut Criterion) {
    let abc = Alphabet::rna();
    let mut cm = load_cm();
    configure_scores(&mut cm); // glocal
    let emap = EmitMap::build(&cm);
    let dsq = read_seq(&abc, "trna_embedded.fa"); // 181 nt, real tRNA at 61..131
    let w = cm.w as usize;

    let mut g = c.benchmark_group("cm_dp");
    g.bench_function("cyk_scan", |b| {
        b.iter(|| cyk_scan(black_box(&cm), black_box(&dsq), w))
    });
    g.bench_function("inside_scan", |b| {
        b.iter(|| inside_scan(black_box(&cm), black_box(&dsq), w))
    });
    g.bench_function("align_glocal", |b| {
        b.iter(|| align_glocal(black_box(&cm), black_box(&dsq), 61, 131, black_box(&emap)))
    });
    g.finish();
}

fn bench_cm_search(c: &mut Criterion) {
    let abc = Alphabet::rna();
    let consensus = {
        let recs = read_fasta(data("trna_cons.fa")).unwrap();
        recs[0].seq.clone()
    };
    let dsq = synthetic_with_trnas(&abc, &consensus, 400);

    let mut cm = load_cm();
    configure_local(&mut cm); // search uses local mode (cmsearch default)
    let w = cm.w as usize;

    let mut g = c.benchmark_group("cm_search");
    g.sample_size(10); // the full-DP CYK scan is heavy (O(L·W²·M)); few samples
    g.bench_function("cyk_search_unfiltered", |b| {
        b.iter(|| cyk_search(black_box(&cm), black_box(&dsq), w, 0.0))
    });
    g.bench_function("cm_pipeline_filtered", |b| {
        b.iter(|| {
            cm_pipeline_search(
                black_box(&cm),
                black_box(&dsq),
                w,
                0.0,
                PipelineParams::default(),
            )
        })
    });
    g.finish();
}

criterion_group!(benches, bench_cm_dp, bench_cm_search);
criterion_main!(benches);
