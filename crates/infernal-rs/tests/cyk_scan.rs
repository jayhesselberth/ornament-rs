//! Differential test: the native glocal CYK scan must reproduce the `cmsearch` oracle's
//! bit score and hit coordinates on the calibrated RF00005 tRNA model.
//!
//! Oracle (Infernal in the pixi `dev` env):
//!   cmsearch -g --max --cyk --toponly --tblout - tRNA.cm trna_cons.fa
//! reports the 71-nt consensus as a hit over seq 1..71 at 87.3 bits (glocal CYK, no
//! filters/bands). The native scanner targets the same number.

use easel_rs::{read_fasta, Alphabet};
use infernal_rs::search::{cyk_scan_glocal, inside_scan_glocal};
use infernal_rs::{configure_scores, evalue, parse_cm_file, SearchMode};

fn cm() -> infernal_rs::Cm {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
    configure_scores(&mut cm);
    cm
}

#[test]
fn glocal_cyk_matches_oracle_on_consensus() {
    let cm = cm();
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_cons.fa");
    let recs = read_fasta(fa).expect("read consensus fasta");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");
    let l = recs[0].seq.len();
    assert_eq!(l, 71, "consensus length");

    let scan = cyk_scan_glocal(&cm, &dsq, cm.w as usize);
    let best = scan.best.expect("a hit");

    eprintln!("native best: score={:.3} i={} j={}", best.score, best.i, best.j);
    // Oracle: 87.3 bits over seq 1..71.
    assert!(
        (best.score - 87.3).abs() < 0.2,
        "native CYK score {:.3} != oracle 87.3",
        best.score
    );
    assert_eq!((best.i, best.j), (1, 71), "hit coordinates");

    // E-value parity: oracle reports 1.1e-17 (glocal CYK, top-only => L searched residues).
    let e = evalue(&cm, SearchMode::GlocalCyk, best.score, l as f64).unwrap();
    eprintln!("native CYK E-value: {e:.2e}");
    assert!(
        (e.log10() - 1.1e-17_f64.log10()).abs() < 0.15,
        "native E-value {e:.3e} != oracle 1.1e-17"
    );
}

#[test]
fn glocal_cyk_finds_embedded_trna() {
    // The tRNA consensus embedded in random flanks at seq 61..131. The scanner must locate
    // it (correct endpoints) and reproduce the oracle's 87.3-bit score, exercising the
    // per-end-position scan rather than a whole-sequence alignment.
    let cm = cm();
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_embedded.fa");
    let recs = read_fasta(fa).expect("read embedded fasta");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");

    let scan = cyk_scan_glocal(&cm, &dsq, cm.w as usize);
    let best = scan.best.expect("a hit");
    eprintln!("embedded best: score={:.3} i={} j={}", best.score, best.i, best.j);

    // Oracle: seq 61..131 at 87.3 bits.
    assert_eq!((best.i, best.j), (61, 131), "located the embedded tRNA");
    assert!(
        (best.score - 87.3).abs() < 0.2,
        "native CYK score {:.3} != oracle 87.3",
        best.score
    );
}

#[test]
fn glocal_inside_matches_oracle() {
    // Oracle (Inside, no --cyk): cmsearch -g --max reports 87.3 bits for both the
    // consensus and the embedded tRNA. Inside >= CYK; here the dominant parse makes them
    // round to the same value. Validates the log-sum-exp semiring.
    let cm = cm();
    let abc = Alphabet::rna();
    for (file, coords) in [
        ("/tests/data/trna_cons.fa", (1usize, 71usize)),
        ("/tests/data/trna_embedded.fa", (61, 131)),
    ] {
        let path = format!("{}{}", env!("CARGO_MANIFEST_DIR"), file);
        let recs = read_fasta(&path).expect("read fasta");
        let dsq = abc.digitize(&recs[0].seq).expect("digitize");
        let scan = inside_scan_glocal(&cm, &dsq, cm.w as usize);
        let best = scan.best.expect("a hit");
        eprintln!("inside {file}: score={:.3} i={} j={}", best.score, best.i, best.j);
        assert_eq!((best.i, best.j), coords, "{file} coordinates");
        assert!(
            (best.score - 87.3).abs() < 0.25,
            "{file} native Inside {:.3} != oracle 87.3",
            best.score
        );
        // Inside must be >= CYK for the same endpoint.
        let cyk = cyk_scan_glocal(&cm, &dsq, cm.w as usize).best.unwrap();
        assert!(best.score >= cyk.score - 1e-3, "Inside < CYK");
    }
}
