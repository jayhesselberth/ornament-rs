//! The p7 Forward-filter windowing pipeline must not drop true hits relative to the
//! brute-force whole-strand scan, while restricting the CM DP to filter survivors.
//!
//! Two invariants are checked:
//!  1. **Parity** — on a sequence with strong tRNAs, `cm_pipeline_search` returns exactly the
//!     same above-cutoff hits (coords/strand/score/E-value) as the unfiltered `cyk_search`
//!     (the `cmsearch --max` equivalent already validated in `multihit_strand.rs`).
//!  2. **Pruning** — on a long sequence with sparse signal, the filter forwards only a small
//!     fraction of residues to the CM, demonstrating the intended speed-up.

use easel_rs::{read_fasta, Alphabet, Dsq};
use infernal_rs::{
    cm_pipeline_search, configure_local, cyk_search, parse_cm_file, Hit, PipelineParams, Strand,
};

fn load_cm() -> infernal_rs::Cm {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
    configure_local(&mut cm); // cmsearch default (local, filters)
    cm
}

/// Match two hit sets ignoring order: same coordinates/strand, scores within Infernal slack.
fn assert_same_hits(got: &[Hit], want: &[Hit]) {
    assert_eq!(
        got.len(),
        want.len(),
        "hit count: pipeline {} vs brute {}",
        got.len(),
        want.len()
    );
    for h in want {
        let m = got
            .iter()
            .find(|g| g.i == h.i && g.j == h.j && g.strand == h.strand)
            .unwrap_or_else(|| panic!("pipeline missing hit {:?} {}..{}", h.strand, h.i, h.j));
        assert!(
            (m.score - h.score).abs() < 1e-3,
            "score at {}..{}: pipeline {} vs brute {}",
            h.i,
            h.j,
            m.score,
            h.score
        );
        match (m.evalue, h.evalue) {
            (Some(a), Some(b)) => assert!(
                (a.log10() - b.log10()).abs() < 1e-6,
                "E-value at {}..{}: {a:.3e} vs {b:.3e}",
                h.i,
                h.j
            ),
            (a, b) => assert_eq!(a.is_some(), b.is_some(), "E-value presence mismatch"),
        }
    }
}

#[test]
fn pipeline_matches_bruteforce_on_multi() {
    let cm = load_cm();
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_multi.fa");
    let recs = read_fasta(fa).expect("read multi fasta");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");

    let cutoff = 20.0;
    let brute = cyk_search(&cm, &dsq, cm.w as usize, cutoff);
    let (filtered, stats) =
        cm_pipeline_search(&cm, &dsq, cm.w as usize, cutoff, PipelineParams::default())
            .expect("pipeline");

    for h in &filtered {
        eprintln!(
            "pipeline: {:?} {}..{} score={:.2} E={:?}",
            h.strand, h.i, h.j, h.score, h.evalue
        );
    }
    eprintln!(
        "windows={} survivors={} residues_to_cm={}/{} ({:.1}%)",
        stats.n_windows,
        stats.n_survivors,
        stats.residues_to_cm,
        stats.residues_searched,
        100.0 * stats.cm_fraction()
    );

    // The three strong tRNAs (2 plus, 1 minus) must survive, bit-identical to the brute scan.
    assert_eq!(filtered.len(), 3, "three significant hits preserved");
    assert_same_hits(&filtered, &brute);
}

/// Deterministic LCG → a random RNA flank that won't spuriously trip the filter at F3=1e-5.
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

/// Build the RF00005 consensus and its reverse complement embedded in long random flanks,
/// so both strands carry one strong (~87-bit) hit, far apart, amid signal-free sequence.
fn synthetic_sparse(abc: &Alphabet, flank: usize) -> (String, String) {
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
    let seq = format!(
        "{}{}{}{}{}",
        random_flank(flank, &mut seed),
        cons,
        random_flank(flank, &mut seed),
        rc_cons,
        random_flank(flank, &mut seed),
    );
    (seq, cons)
}

/// Pruning + no-drop *with the filter active* (multiple tiles): on a moderate sparse
/// sequence the pipeline matches the brute-force scan exactly while skipping CM work on the
/// signal-free flanks.
#[test]
fn pipeline_prunes_and_preserves() {
    let cm = load_cm();
    let abc = Alphabet::rna();
    let (seq, _) = synthetic_sparse(&abc, 300);
    let dsq = abc.digitize(&seq).expect("digitize");

    let cutoff = 25.0;
    let brute = cyk_search(&cm, &dsq, cm.w as usize, cutoff);
    let (filtered, stats) =
        cm_pipeline_search(&cm, &dsq, cm.w as usize, cutoff, PipelineParams::default())
            .expect("pipeline");

    eprintln!(
        "L={} windows={} survivors={} residues_to_cm={}/{} ({:.1}%)",
        seq.len(),
        stats.n_windows,
        stats.n_survivors,
        stats.residues_to_cm,
        stats.residues_searched,
        100.0 * stats.cm_fraction()
    );

    assert_same_hits(&filtered, &brute); // no true hit dropped, scores/coords identical
    assert!(
        filtered.iter().any(|h| h.strand == Strand::Plus)
            && filtered.iter().any(|h| h.strand == Strand::Minus),
        "one strong hit on each strand"
    );
    // The filter genuinely pruned residues (some flank tiles failed F3).
    assert!(
        stats.residues_to_cm < stats.residues_searched,
        "expected pruning, sent all {} residues to CM",
        stats.residues_to_cm
    );
}

/// Scale demonstration: on a long sparse sequence the pipeline forwards only a small fraction
/// of residues to the (expensive) CM, yet still recovers both embedded tRNAs. Runs the
/// pipeline alone — the brute-force scan over the full `2L` is what we are accelerating, and
/// is far too slow here; its no-drop equivalence is established by the tests above.
#[test]
fn pipeline_scales_on_long_sparse_sequence() {
    let cm = load_cm();
    let abc = Alphabet::rna();
    let (seq, _) = synthetic_sparse(&abc, 2500); // ~7.6 kb, two tRNAs
    let dsq = abc.digitize(&seq).expect("digitize");

    let cutoff = 25.0;
    let (filtered, stats) =
        cm_pipeline_search(&cm, &dsq, cm.w as usize, cutoff, PipelineParams::default())
            .expect("pipeline");

    eprintln!(
        "L={} windows={} survivors={} residues_to_cm={}/{} ({:.1}%)",
        seq.len(),
        stats.n_windows,
        stats.n_survivors,
        stats.residues_to_cm,
        stats.residues_searched,
        100.0 * stats.cm_fraction()
    );
    for h in &filtered {
        eprintln!("hit: {:?} {}..{} score={:.2}", h.strand, h.i, h.j, h.score);
    }

    // Exactly the two embedded tRNAs (one per strand), each at full CM score.
    assert_eq!(filtered.len(), 2, "the two embedded tRNAs");
    assert!(filtered.iter().all(|h| h.score > 70.0), "both strong hits");
    assert_eq!(
        filtered.iter().filter(|h| h.strand == Strand::Plus).count(),
        1
    );
    assert_eq!(
        filtered
            .iter()
            .filter(|h| h.strand == Strand::Minus)
            .count(),
        1
    );
    // The CM saw only a small fraction of the searched residues.
    assert!(
        stats.cm_fraction() < 0.25,
        "filter should prune most residues, got {:.1}%",
        100.0 * stats.cm_fraction()
    );
}
