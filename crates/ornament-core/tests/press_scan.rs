//! End-to-end test that a pressed CM database produces identical scan hits to the `.cm` it was
//! built from — the guarantee the `press` optimization rests on.

use ornament_core::infernal::{press_cm, scan_native_multi, CMHit};

const CM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
const FA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_multi.fa");
const E_VALUE: f64 = 1e-5;

/// Comparable projection of a hit (E-value/score are floats; compare coordinates + strand + names).
fn key(h: &CMHit) -> (String, String, usize, usize, char) {
    (
        h.query_name.clone(),
        h.target_name.clone(),
        h.target_start,
        h.target_end,
        h.strand,
    )
}

#[test]
fn pressed_scan_matches_unpressed_scan() {
    let tmp = tempfile::tempdir().unwrap();
    let orm = tmp.path().join("tRNA.cm.orm");

    let (n, written) = press_cm(CM, Some(orm.as_path()), false).unwrap();
    assert_eq!(n, 1);
    assert_eq!(written, orm);

    let mut from_cm = scan_native_multi(CM, FA, E_VALUE).unwrap();
    let mut from_orm = scan_native_multi(orm.to_str().unwrap(), FA, E_VALUE).unwrap();
    assert!(!from_cm.is_empty(), "fixture should produce hits");

    // Hits are already best-first sorted by both calls; compare the coordinate/strand projection.
    let a: Vec<_> = from_cm.drain(..).map(|h| key(&h)).collect();
    let b: Vec<_> = from_orm.drain(..).map(|h| key(&h)).collect();
    assert_eq!(
        a, b,
        "pressed scan must return identical hits to the .cm scan"
    );
}
