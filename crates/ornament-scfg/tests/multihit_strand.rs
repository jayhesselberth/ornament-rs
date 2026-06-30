//! Differential test for both-strand, multi-hit scanning with gamma overlap resolution.
//!
//! Oracle (Infernal in the pixi `dev` env), LOCAL mode (cmsearch default), no filters:
//!   cmsearch --max --cyk --tblout trna_multi.cyk.tblout tRNA.cm trna_multi.fa
//! The fixture `trna_multi.fa` embeds the RF00005 consensus tRNA three times in random
//! flanks: two on the plus strand and one reverse-complemented (minus strand). The oracle
//! reports all three at 87.2 bits (E=5e-30), the minus-strand copy as `seq from > seq to`,
//! plus several sub-cutoff spurious hits. The native `cyk_search` must reproduce the three
//! significant hits on the correct strands, coordinates, scores, and E-values.

use ornament_alphabet::{read_fasta, Alphabet};
use ornament_scfg::{
    calc_qdb_bands, configure_local, cyk_search, cyk_search_banded, parse_cm_file, QdbBands, Strand,
};

#[derive(Debug)]
struct Oracle {
    from: usize,
    to: usize,
    strand: char,
    score: f32,
    evalue: f64,
}

fn read_tblout(path: &str) -> Vec<Oracle> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            // cols: ... seqfrom(7) seqto(8) strand(9) ... score(14) evalue(15)
            Oracle {
                from: f[7].parse().unwrap(),
                to: f[8].parse().unwrap(),
                strand: f[9].chars().next().unwrap(),
                score: f[14].parse().unwrap(),
                evalue: f[15].parse().unwrap(),
            }
        })
        .collect()
}

#[test]
fn multihit_both_strands_match_oracle() {
    let cm_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(cm_path).expect("parse tRNA.cm");
    configure_local(&mut cm); // cmsearch default

    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_multi.fa");
    let recs = read_fasta(fa).expect("read multi fasta");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");

    // Report hits >= 20 bits: keeps the three strong tRNAs, drops the spurious < 12-bit hits.
    let hits = cyk_search(&cm, &dsq, cm.w as usize, 20.0);
    for h in &hits {
        eprintln!(
            "native: {:?} {}..{} score={:.2} E={:?}",
            h.strand, h.i, h.j, h.score, h.evalue
        );
    }
    assert_eq!(hits.len(), 3, "three significant hits");

    // Oracle significant hits (score >= 20).
    let oracle: Vec<Oracle> = read_tblout(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/trna_multi.cyk.tblout"
    ))
    .into_iter()
    .filter(|o| o.score >= 20.0)
    .collect();
    assert_eq!(oracle.len(), 3, "three significant oracle hits");

    // Match each native hit to an oracle hit by coordinates (i==from, j==to).
    for h in &hits {
        let nstrand = if h.strand == Strand::Plus { '+' } else { '-' };
        let m = oracle
            .iter()
            .find(|o| o.from == h.i && o.to == h.j)
            .unwrap_or_else(|| panic!("no oracle hit at {}..{}", h.i, h.j));
        assert_eq!(nstrand, m.strand, "strand at {}..{}", h.i, h.j);
        assert!(
            (h.score - m.score).abs() < 0.2,
            "score {} vs oracle {} at {}..{}",
            h.score,
            m.score,
            h.i,
            h.j
        );
        let e = h.evalue.expect("calibrated E-value");
        assert!(
            (e.log10() - m.evalue.log10()).abs() < 0.2,
            "E-value {:.3e} vs oracle {:.3e} at {}..{}",
            e,
            m.evalue,
            h.i,
            h.j
        );
    }

    // One hit on the minus strand (the reverse-complemented copy), reported high->low.
    let minus: Vec<_> = hits.iter().filter(|h| h.strand == Strand::Minus).collect();
    assert_eq!(minus.len(), 1, "exactly one minus-strand hit");
    assert!(minus[0].i > minus[0].j, "minus hit reported as from > to");
    // Two hits on the plus strand (multi-hit on one strand).
    assert_eq!(
        hits.iter().filter(|h| h.strand == Strand::Plus).count(),
        2,
        "two plus-strand hits"
    );

    // Sorted best-first by E-value.
    for w in hits.windows(2) {
        assert!(
            w[0].evalue.unwrap() <= w[1].evalue.unwrap(),
            "sorted by E-value"
        );
    }
}

/// The query-dependent-banded brute search must return exactly what the unbanded brute search
/// returns — the bands keep every real hit's optimal parse in range. This is the non-pipeline
/// counterpart of `pipeline::pipeline_matches_bruteforce` (which proves it for the windowing
/// path), and it guards `cyk_search_banded` (used by `ornament-core`'s native scan).
#[test]
fn banded_brute_matches_unbanded() {
    let cm_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(cm_path).expect("parse tRNA.cm");
    configure_local(&mut cm);

    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_multi.fa");
    let dsq = abc
        .digitize(&read_fasta(fa).expect("read fasta")[0].seq)
        .expect("digitize");

    let w = cm.w as usize;
    let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);
    let unbanded = cyk_search(&cm, &dsq, w, 20.0);
    let banded = cyk_search_banded(&cm, &dsq, w, 20.0, &bands);

    assert!(!unbanded.is_empty(), "expected hits to compare");
    assert_eq!(
        banded, unbanded,
        "banded brute search differs from unbanded"
    );
}
