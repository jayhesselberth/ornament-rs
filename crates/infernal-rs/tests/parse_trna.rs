//! Integration test: parse the calibrated RF00005 tRNA `.cm` fixture and check parity
//! with its header values and a couple of hand-verified state lines.

use infernal_rs::model::{exp_mode, st, stid};
use infernal_rs::parse_cm_file;

fn fixture() -> infernal_rs::Cm {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    parse_cm_file(path).expect("parse tRNA.cm")
}

#[test]
fn header_fields_match() {
    let cm = fixture();
    assert_eq!(cm.name, "tRNA");
    assert_eq!(cm.acc.as_deref(), Some("RF00005"));
    assert_eq!(cm.m, 227, "STATES");
    assert_eq!(cm.nodes, 60, "NODES");
    assert_eq!(cm.clen, 71, "CLEN (recomputed from MATP/MATL/MATR nodes)");
    assert_eq!(cm.w, 218, "W");
    assert_eq!(cm.k(), 4);
    // All 227 state lines parsed and validated; arrays sized to M.
    assert_eq!(cm.sttype.len(), 227);
    assert_eq!(cm.t.len(), 227);
}

#[test]
fn evalue_exp_tails_parsed() {
    let cm = fixture();
    // ECMLC    0.79302    -4.21089     2.90112     1600000      337752  0.003553
    let lc = cm.exp[exp_mode::CM_LC];
    assert!(lc.is_valid);
    assert!((lc.lambda - 0.79302).abs() < 1e-5);
    assert!((lc.mu_orig - 2.90112).abs() < 1e-5);
    assert_eq!(lc.nrandhits, 337752);
    assert!((lc.dbsize - 1_600_000.0).abs() < 1.0);
    // ECMGI    0.41538   -11.77572    -1.10020     1600000       33721  0.011862
    let gi = cm.exp[exp_mode::CM_GI];
    assert!((gi.lambda - 0.41538).abs() < 1e-5);
    assert_eq!(gi.nrandhits, 33721);
}

#[test]
fn first_states_have_expected_topology() {
    let cm = fixture();
    // State 0: ROOT_S start, 4 children, no emissions.
    assert_eq!(cm.sttype[0], st::S);
    assert_eq!(cm.stid[0], stid::ROOT_S);
    assert_eq!(cm.cnum[0], 4);
    assert_eq!(cm.cfirst[0], 1);
    assert_eq!(cm.t[0].len(), 4);
    assert!(cm.e[0].is_empty());

    // State 6 is the first MATP MP (16 pair emissions, 6 transitions).
    assert_eq!(cm.sttype[6], st::MP);
    assert_eq!(cm.stid[6], stid::MATP_MP);
    assert_eq!(cm.e[6].len(), 16);
    assert_eq!(cm.t[6].len(), cm.cnum[6] as usize);
    // Probabilities are normalized: emissions sum to ~1, transitions sum to ~1.
    let esum: f32 = cm.e[6].iter().sum();
    assert!((esum - 1.0).abs() < 1e-2, "MP emissions sum = {esum}");
    let tsum: f32 = cm.t[6].iter().sum();
    assert!((tsum - 1.0).abs() < 1e-2, "MP transitions sum = {tsum}");
}

#[test]
fn configure_scores_round_trips_file_values() {
    // The .cm stores log-odds scores; parse converts to probabilities; configure converts
    // back. The round trip must reproduce the on-disk scores to float precision.
    let mut cm = fixture();
    infernal_rs::configure_scores(&mut cm);

    // State 0 ROOT_S transitions on disk: -10.978 -12.223 -0.031 -5.613
    let want_t0 = [-10.978f32, -12.223, -0.031, -5.613];
    for (i, w) in want_t0.iter().enumerate() {
        assert!(
            (cm.tsc[0][i] - w).abs() < 1e-2,
            "tsc[0][{i}] = {} want {w}",
            cm.tsc[0][i]
        );
    }

    // State 6 MATP_MP pair emission scores on disk begin: -5.257 -3.238 -3.838 1.469 ...
    let want_e6 = [-5.257f32, -3.238, -3.838, 1.469];
    for (i, w) in want_e6.iter().enumerate() {
        assert!(
            (cm.esc[6][i] - w).abs() < 1e-2,
            "esc[6][{i}] = {} want {w}",
            cm.esc[6][i]
        );
    }
}

#[test]
fn nodemap_and_ndidx_are_consistent() {
    let cm = fixture();
    // Node 0 is ROOT, its first state is 0.
    assert_eq!(cm.nodemap[0], 0);
    assert_eq!(cm.ndidx[0], 0);
    // Every state's node index is in range, and nodemap points back consistently.
    for v in 0..cm.m {
        let nd = cm.ndidx[v];
        assert!(nd < cm.nodes);
        assert!(cm.nodemap[nd] <= v);
    }
    // Embedded HMMER3 filter HMM was captured.
    assert!(cm.fp7_text.as_deref().unwrap_or("").contains("HMMER3/f"));
}
