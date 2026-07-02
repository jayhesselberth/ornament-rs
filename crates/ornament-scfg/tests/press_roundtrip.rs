//! Round-trip tests for the pressed-CM database (`press` module).
//!
//! A pressed model must reproduce, bit-for-bit, the state a fresh parse + configure + QDB band
//! computation would produce — that is what lets `scan`/`search` skip the work and still return
//! identical hits. These tests also cover the format guards (magic / version).

use ornament_scfg::{
    calc_qdb_bands, configure_local, configure_scores, load_pressed, parse_cm_records_file,
    press_cm_file, QdbBands,
};

const CM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");

#[test]
fn pressed_model_matches_a_fresh_prep() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("tRNA.cm.orm");

    // Fresh prep of model 0, exactly as the scan pipeline builds it.
    let mut models = parse_cm_records_file(CM).unwrap();
    let mut cm = models.remove(0);
    let mut align_cm = cm.clone();
    configure_scores(&mut align_cm);
    configure_local(&mut cm);
    let scan_bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);
    let align_bands = calc_qdb_bands(&align_cm, QdbBands::DEFAULT_BETA);

    // Press + reload.
    let n = press_cm_file(std::path::Path::new(CM), &out, false).unwrap();
    assert_eq!(n, 1, "fixture has one model");
    let db = load_pressed(&out).unwrap();
    assert_eq!(db.models.len(), 1);
    let pm = &db.models[0];

    // Identity metadata + topology survive the round trip.
    assert_eq!(pm.cm.name, cm.name);
    assert_eq!(pm.cm.m, cm.m);
    assert_eq!(pm.cm.clen, cm.clen);
    assert_eq!(pm.w_max, cm.w as usize);

    // The alphabet was reattached (not the deserialization placeholder) and matches.
    assert_eq!(pm.cm.abc.kp, cm.abc.kp);
    assert_eq!(pm.cm.abc.k, cm.abc.k);
    assert!(
        std::sync::Arc::ptr_eq(&pm.cm.abc, &pm.align_cm.abc),
        "one shared alphabet"
    );

    // Log-odds scores match (config was captured, not recomputed differently).
    assert_eq!(pm.cm.tsc, cm.tsc, "local scan scores round-trip exactly");
    assert_eq!(
        pm.align_cm.tsc, align_cm.tsc,
        "glocal align scores round-trip exactly"
    );

    // The stored bands equal a fresh computation on the corresponding models.
    assert_eq!(pm.scan_bands.dmin, scan_bands.dmin);
    assert_eq!(pm.scan_bands.dmax, scan_bands.dmax);
    assert_eq!(pm.scan_bands.w, scan_bands.w);
    assert_eq!(pm.align_bands.dmin, align_bands.dmin);
    assert_eq!(pm.align_bands.dmax, align_bands.dmax);

    // The filter cascade rebuilds from the stored HMM + bands without error.
    assert!(pm.prepared_filters().is_ok());
}

#[test]
fn load_rejects_non_pressed_file() {
    // The plain `.cm` text is not a pressed DB — the loader must reject it, not decode garbage.
    let err = load_pressed(std::path::Path::new(CM)).unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("not a valid pressed DB") || msg.contains("not a pressed"),
        "unexpected error: {msg}"
    );
}

#[test]
fn press_refuses_to_overwrite_without_force() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("tRNA.cm.orm");
    press_cm_file(std::path::Path::new(CM), &out, false).unwrap();
    // Second press without force must fail; with force it succeeds.
    assert!(press_cm_file(std::path::Path::new(CM), &out, false).is_err());
    assert!(press_cm_file(std::path::Path::new(CM), &out, true).is_ok());
}
