//! The rayon-parallel search must be *deterministic*: its output is independent of how many
//! threads run it. A single-thread rayon pool is the scalar reference (every window/strand
//! computed sequentially); a many-thread pool exercises the fork-join and the parallel tile
//! sweep + per-survivor-window scans. Both must return bit-identical hits.
//!
//! This is the "SIMD/rayon output == scalar reference" invariant from issue #5 for the rayon
//! layer: the per-window DP is unchanged, and results are merged in a fixed order then sorted,
//! so thread count must never perturb scores, coordinates, strand, or E-values.

use easel_rs::{read_fasta, Alphabet, Dsq};
use infernal_rs::{
    cm_pipeline_search, configure_local, cyk_search, parse_cm_file, Hit, PipelineParams,
};

fn load_cm() -> infernal_rs::Cm {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
    configure_local(&mut cm);
    cm
}

/// Deterministic LCG → random RNA flank (mirrors `pipeline.rs`).
fn random_flank(n: usize, seed: &mut u64) -> String {
    let bases = [b'A', b'C', b'G', b'U'];
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((*seed >> 33) & 3) as usize] as char
        })
        .collect()
}

/// RF00005 consensus + its reverse complement embedded in long random flanks: one strong hit
/// per strand, far apart, spanning many filter tiles so the parallel sweep does real work.
fn synthetic_sparse(abc: &Alphabet, flank: usize) -> String {
    let cons = {
        let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_cons.fa");
        read_fasta(fa).unwrap()[0].seq.clone()
    };
    let rc_cons = {
        let mut d: Vec<Dsq> = abc.digitize(&cons).unwrap();
        abc.revcomp(&mut d).unwrap();
        abc.textize(&d)
    };
    let mut seed = 0x0123_4567_89ab_cdefu64;
    format!(
        "{}{}{}{}{}",
        random_flank(flank, &mut seed),
        cons,
        random_flank(flank, &mut seed),
        rc_cons,
        random_flank(flank, &mut seed),
    )
}

/// Run `f` inside a rayon pool with exactly `threads` worker threads.
fn with_threads<T: Send>(threads: usize, f: impl FnOnce() -> T + Send) -> T {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build rayon pool");
    pool.install(f)
}

#[test]
fn cyk_search_is_thread_count_independent() {
    let cm = load_cm();
    let abc = Alphabet::rna();
    let dsq = abc
        .digitize(&synthetic_sparse(&abc, 200))
        .expect("digitize");
    let w = cm.w as usize;

    let one: Vec<Hit> = with_threads(1, || cyk_search(&cm, &dsq, w, 25.0));
    let many: Vec<Hit> = with_threads(8, || cyk_search(&cm, &dsq, w, 25.0));

    assert!(!one.is_empty(), "expected hits to compare");
    assert_eq!(one, many, "cyk_search differs by thread count");
}

#[test]
fn pipeline_is_thread_count_independent() {
    let cm = load_cm();
    let abc = Alphabet::rna();
    let dsq = abc
        .digitize(&synthetic_sparse(&abc, 400))
        .expect("digitize");
    let w = cm.w as usize;
    let p = PipelineParams::default();

    let (one_hits, one_stats) =
        with_threads(1, || cm_pipeline_search(&cm, &dsq, w, 25.0, p)).expect("pipeline 1t");
    let (many_hits, many_stats) =
        with_threads(8, || cm_pipeline_search(&cm, &dsq, w, 25.0, p)).expect("pipeline 8t");

    assert!(!one_hits.is_empty(), "expected hits to compare");
    assert_eq!(one_hits, many_hits, "pipeline hits differ by thread count");
    // The pruning bookkeeping must also be identical (tile sweep is order-independent).
    assert_eq!(one_stats.n_windows, many_stats.n_windows);
    assert_eq!(one_stats.n_survivors, many_stats.n_survivors);
    assert_eq!(one_stats.residues_to_cm, many_stats.residues_to_cm);
}
