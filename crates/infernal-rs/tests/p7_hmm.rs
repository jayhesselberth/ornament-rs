//! Parse the HMMER3/f filter HMM embedded in the RF00005 `.cm` and check it against the
//! known header values + probability invariants.

use easel_rs::{read_fasta, Alphabet};
use hmmer_rs::hmm::tr;
use hmmer_rs::{forward_bits, parse_p7_hmm};
use infernal_rs::parse_cm_file;

#[test]
fn parses_embedded_filter_hmm() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let cm = parse_cm_file(path).expect("parse cm");
    let text = cm.fp7_text.as_deref().expect("embedded HMM captured");

    let hmm = parse_p7_hmm(text).expect("parse p7 HMM");
    assert_eq!(hmm.name, "tRNA");
    assert_eq!(hmm.m, 71, "LENG");
    assert_eq!(hmm.k, 4);
    assert_eq!(hmm.maxl, 253);

    // Calibration STATS lines.
    assert!((hmm.msv.mu - -8.9378).abs() < 1e-3);
    assert!((hmm.msv.lambda - 0.73469).abs() < 1e-4);
    assert!((hmm.vit.mu - -9.6821).abs() < 1e-3);
    assert!((hmm.fwd.mu - -2.5854).abs() < 1e-3);

    // Match emissions are probabilities summing to ~1 at every node.
    for k in 1..=hmm.m {
        let s: f32 = hmm.mat[k].iter().sum();
        assert!((s - 1.0).abs() < 1e-2, "node {k} match emissions sum = {s}");
    }
    // Insert emissions at node 1 are uniform (0.25 each: -ln 0.25 = 1.38629).
    for &p in &hmm.ins[1] {
        assert!((p - 0.25).abs() < 1e-3);
    }
    // Node transitions are probabilities; M->{M,I,D} sums to ~1.
    let tsum = hmm.t[1][tr::MM] + hmm.t[1][tr::MI] + hmm.t[1][tr::MD];
    assert!((tsum - 1.0).abs() < 1e-2, "node1 M-transitions sum = {tsum}");

    // COMPO present.
    assert!(hmm.compo.is_some());
}

#[test]
fn forward_filter_score_matches_hmmsearch() {
    // The native p7 Forward filter must reproduce hmmsearch's full-sequence Forward bit
    // score (hmmsearch --max on the extracted filter HMM): 53.5 on the consensus, 50.1 on
    // the embedded tRNA.
    let cm = parse_cm_file(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm")).unwrap();
    let hmm = parse_p7_hmm(cm.fp7_text.as_deref().unwrap()).unwrap();
    let abc = Alphabet::rna();

    for (file, want) in [
        ("/tests/data/trna_cons.fa", 53.5f32),
        ("/tests/data/trna_embedded.fa", 50.1f32),
    ] {
        let path = format!("{}{}", env!("CARGO_MANIFEST_DIR"), file);
        let recs = read_fasta(&path).unwrap();
        let dsq = abc.digitize(&recs[0].seq).unwrap();
        let bits = forward_bits(&hmm, &abc, &dsq);
        eprintln!("{file}: native Forward = {bits:.2} bits (oracle {want})");
        assert!(
            (bits - want).abs() < 0.5,
            "{file}: native Forward {bits:.3} != hmmsearch {want}"
        );
    }
}
