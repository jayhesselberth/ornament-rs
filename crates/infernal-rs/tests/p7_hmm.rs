//! Parse the HMMER3/f filter HMM embedded in the RF00005 `.cm` and check it against the
//! known header values + probability invariants.

use easel_rs::{read_fasta, Alphabet};
use hmmer_rs::hmm::tr;
use hmmer_rs::profile::{bg_freqs, p7p, P7Profile};
use hmmer_rs::{
    forward_bits, forward_nats, msv_bits, msv_nats, parse_p7_hmm, viterbi_bits, viterbi_nats,
};
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
    assert!(
        (tsum - 1.0).abs() < 1e-2,
        "node1 M-transitions sum = {tsum}"
    );

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

/// The one rigorous cross-filter invariant is `Viterbi ≤ Forward`: both run on the *same*
/// configured profile, and Viterbi's max-over-paths can never exceed Forward's
/// sum-over-paths (`max ≤ log-sum-exp`). `hmmsearch` only reports the Forward score, so for
/// Viterbi we lean on this bound plus sane-value checks; MSV gets its own exact oracle in
/// [`msv_equals_viterbi_on_msv_profile`].
///
/// NOTE: the naive intuition "MSV ≤ Viterbi ≤ Forward" does **not** hold, and we
/// deliberately do not assert it. `p7_GMSV` waives every core transition penalty (all
/// `t_MM = 0`) and uses a flat `log(2/(M(M+1)))` entry instead of the occupancy-weighted
/// local entry, so on a well-matching target MSV routinely *outscores* Viterbi and even
/// Forward (it is a different, more permissive model). This is exactly why HMMER calibrates
/// MSV/Viterbi/Forward with separate Gumbel parameters and compares each filter to its own
/// threshold, never to one another — see the `STATS LOCAL {MSV,VITERBI,FORWARD}` lines.
#[test]
fn filter_scores_are_sane() {
    let cm = parse_cm_file(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm")).unwrap();
    let hmm = parse_p7_hmm(cm.fp7_text.as_deref().unwrap()).unwrap();
    let abc = Alphabet::rna();

    let eps = 1e-3f32; // f32 slack on the max ≤ log-sum-exp bound
    for file in ["/tests/data/trna_cons.fa", "/tests/data/trna_embedded.fa"] {
        let path = format!("{}{}", env!("CARGO_MANIFEST_DIR"), file);
        let recs = read_fasta(&path).unwrap();
        let dsq = abc.digitize(&recs[0].seq).unwrap();

        let msv = msv_bits(&hmm, &abc, &dsq);
        let vit = viterbi_bits(&hmm, &abc, &dsq);
        let fwd = forward_bits(&hmm, &abc, &dsq);
        eprintln!("{file}: MSV = {msv:.2}  Viterbi = {vit:.2}  Forward = {fwd:.2} bits");

        assert!(msv.is_finite() && vit.is_finite() && fwd.is_finite());
        // The guaranteed cross-filter bound (same profile, max ≤ sum).
        assert!(vit <= fwd + eps, "{file}: Viterbi {vit:.3} > Forward {fwd:.3}");
        // All three should land well above the null on a true tRNA.
        assert!(
            msv > 0.0 && vit > 0.0 && fwd > 0.0,
            "{file}: expected positive filter scores, got MSV {msv:.3} Vit {vit:.3} Fwd {fwd:.3}"
        );
    }
}

/// Exact internal oracle for MSV, mirroring HMMER's `utest_msv` (`generic_msv.c`): the MSV
/// score equals a *Viterbi* score on a profile whose core transitions are rigged so that
/// only ungapped diagonals survive — all `t_MM = 0`, every other core transition `-inf`,
/// and `t_BMk` uniform at `log(2/(M(M+1)))`. Under our multihit-local config (`nj = 1`,
/// E move/loop `= -log 2`) this rigged-Viterbi reproduces `p7_GMSV` with `nu = 2.0`
/// (their length models coincide: `tmove = log(3/(L+3))`, `tloop = log(L/(L+3))`), so the
/// two raw nat scores must agree to within f32 noise.
#[test]
fn msv_equals_viterbi_on_msv_profile() {
    let cm = parse_cm_file(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm")).unwrap();
    let hmm = parse_p7_hmm(cm.fp7_text.as_deref().unwrap()).unwrap();
    let abc = Alphabet::rna();
    let bg = bg_freqs(hmm.k);

    for file in ["/tests/data/trna_cons.fa", "/tests/data/trna_embedded.fa"] {
        let path = format!("{}{}", env!("CARGO_MANIFEST_DIR"), file);
        let recs = read_fasta(&path).unwrap();
        let dsq = abc.digitize(&recs[0].seq).unwrap();
        let l = dsq.len().saturating_sub(2);

        let prof = P7Profile::config_local(&hmm, &abc, &bg, l);
        let msv = msv_nats(&prof, &dsq);

        // Rig a clone for the "MSV via Viterbi" equivalence (cf. utest_msv).
        let mut g2 = prof.clone();
        let m = g2.m;
        let tbmk = (2.0f32 / (m as f32 * (m as f32 + 1.0))).ln();
        for v in g2.tsc.iter_mut() {
            *v = f32::NEG_INFINITY;
        }
        for k in 1..m {
            g2.tsc[k * p7p::NTRANS + p7p::MM] = 0.0; // M_k -> M_{k+1}
        }
        for k in 0..m {
            g2.tsc[k * p7p::NTRANS + p7p::BM] = tbmk; // B -> M_{k+1}
        }
        let vit = viterbi_nats(&g2, &dsq);

        eprintln!("{file}: MSV nats {msv:.5} vs rigged-Viterbi nats {vit:.5}");
        assert!(
            (msv - vit).abs() < 1e-2,
            "{file}: MSV {msv:.5} != rigged Viterbi {vit:.5}"
        );
    }
}

/// `_nats`/`_bits` wrappers must be consistent: bits = (nats - null_one(L)) / ln2, and the
/// public bit wrappers must agree with re-deriving from the nat scores via a fresh profile.
#[test]
fn nat_and_bit_wrappers_agree() {
    let cm = parse_cm_file(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm")).unwrap();
    let hmm = parse_p7_hmm(cm.fp7_text.as_deref().unwrap()).unwrap();
    let abc = Alphabet::rna();
    let bg = bg_freqs(hmm.k);

    let path = format!("{}{}", env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_cons.fa");
    let recs = read_fasta(&path).unwrap();
    let dsq = abc.digitize(&recs[0].seq).unwrap();
    let l = dsq.len().saturating_sub(2);
    let prof = P7Profile::config_local(&hmm, &abc, &bg, l);

    let null = hmmer_rs::profile::null_one(l);
    let to_bits = |nats: f32| (nats - null) / std::f32::consts::LN_2;

    assert!((to_bits(msv_nats(&prof, &dsq)) - msv_bits(&hmm, &abc, &dsq)).abs() < 1e-3);
    assert!((to_bits(viterbi_nats(&prof, &dsq)) - viterbi_bits(&hmm, &abc, &dsq)).abs() < 1e-3);
    assert!((to_bits(forward_nats(&prof, &dsq)) - forward_bits(&hmm, &abc, &dsq)).abs() < 1e-3);
}
