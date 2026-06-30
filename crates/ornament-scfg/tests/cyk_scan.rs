//! Differential test: the native glocal CYK scan must reproduce the `cmsearch` oracle's
//! bit score and hit coordinates on the calibrated RF00005 tRNA model.
//!
//! Oracle (Infernal in the pixi `dev` env):
//!   cmsearch -g --max --cyk --toponly --tblout - tRNA.cm trna_cons.fa
//! reports the 71-nt consensus as a hit over seq 1..71 at 87.3 bits (glocal CYK, no
//! filters/bands). The native scanner targets the same number.

use ornament_alphabet::{read_fasta, Alphabet};
use ornament_scfg::search::{cyk_scan, cyk_scan_glocal, inside_scan, inside_scan_glocal};
use ornament_scfg::{configure_local, configure_scores, evalue, parse_cm_file, SearchMode};

fn cm() -> ornament_scfg::Cm {
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

    eprintln!(
        "native best: score={:.3} i={} j={}",
        best.score, best.i, best.j
    );
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
    eprintln!(
        "embedded best: score={:.3} i={} j={}",
        best.score, best.i, best.j
    );

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
        eprintln!(
            "inside {file}: score={:.3} i={} j={}",
            best.score, best.i, best.j
        );
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

#[test]
fn local_cyk_and_inside_match_oracle() {
    // Local mode (cmsearch default, no -g). Oracle with --max [--cyk]: the consensus and the
    // embedded tRNA both score 87.2 bits (slightly below glocal's 87.3 — the local begin/end
    // penalty), best hit at the same coordinates.
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse");
    configure_local(&mut cm); // local instead of glocal
    assert!(cm.is_local);
    let abc = Alphabet::rna();

    for (file, coords) in [
        ("/tests/data/trna_cons.fa", (1usize, 71usize)),
        ("/tests/data/trna_embedded.fa", (61, 131)),
    ] {
        let p = format!("{}{}", env!("CARGO_MANIFEST_DIR"), file);
        let recs = read_fasta(&p).expect("read");
        let dsq = abc.digitize(&recs[0].seq).expect("digitize");

        let cyk = cyk_scan(&cm, &dsq, cm.w as usize).best.expect("cyk hit");
        eprintln!("local CYK {file}: {:.3} {}..{}", cyk.score, cyk.i, cyk.j);
        assert_eq!((cyk.i, cyk.j), coords, "{file} local CYK coords");
        assert!(
            (cyk.score - 87.2).abs() < 0.2,
            "{file} local CYK {:.3} != 87.2",
            cyk.score
        );

        let ins = inside_scan(&cm, &dsq, cm.w as usize)
            .best
            .expect("inside hit");
        assert!(
            (ins.score - 87.2).abs() < 0.3,
            "{file} local Inside {:.3} != 87.2",
            ins.score
        );
        assert!(ins.score >= cyk.score - 1e-3);
    }
}

#[test]
fn glocal_config_unaffected_by_local_paths() {
    // A glocal-configured CM must produce exactly the glocal result (local code paths are
    // no-ops when is_local is false / beginsc/endsc are IMPOSSIBLE).
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse");
    configure_scores(&mut cm);
    assert!(!cm.is_local);
    let abc = Alphabet::rna();
    let recs = read_fasta(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/trna_cons.fa"
    ))
    .unwrap();
    let dsq = abc.digitize(&recs[0].seq).unwrap();
    let best = cyk_scan(&cm, &dsq, cm.w as usize).best.unwrap();
    assert!((best.score - 87.3).abs() < 0.2, "glocal still 87.3");
}
