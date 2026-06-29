//! Fast natural-log log-sum-exp via a precomputed table — a port of HMMER's
//! `p7_FLogsumInit`/`p7_FLogsum` (`logsum.c`).
//!
//! `log(e^a + e^b) = max + log(1 + e^-(max-min))`. The correction term `log(1 + e^-d)`
//! depends only on the non-negative difference `d`, so HMMER tabulates it on a fine grid and
//! replaces the per-cell `exp`+`log` (which dominated the Forward DP — ~45% of a genome scan)
//! with a single table lookup. The table is the exact function HMMER ships, so using it makes
//! the native Forward score *more* faithful to `hmmsearch` (which also uses it), not less.
//!
//! Accuracy: the grid step is `1/SCALE` nats and the max absolute error of the correction is
//! well under `0.001` nats (HMMER's documented bound), negligible against the filter's bit
//! thresholds. Differences `≥ MAXDIFF` contribute `< 2^-23` to the float result and are
//! treated as `max` exactly, as in HMMER.

use std::sync::OnceLock;

/// Grid resolution: table index = `round(d * SCALE)` for difference `d` in nats.
const SCALE: f32 = 1000.0;
/// Table length: covers `d ∈ [0, MAXDIFF)` at `1/SCALE` spacing (`16.0 * 1000`).
const TBL: usize = 16000;
/// Above this difference (nats) the correction underflows in `f32`; return `max` directly.
const MAXDIFF: f32 = 15.7;

static TABLE: OnceLock<Box<[f32]>> = OnceLock::new();

/// The lazily-built correction table; `table()[i] = log(1 + e^(-i/SCALE))`. The first call
/// fills it (≈16k `exp`/`log`), every later call is a cheap pointer return. Hoist it out of
/// the DP's hot loop (call once, index the returned slice) to avoid per-cell `OnceLock` cost.
#[inline]
pub fn table() -> &'static [f32] {
    TABLE.get_or_init(|| {
        (0..TBL)
            .map(|i| (1.0 + (-(i as f32) / SCALE).exp()).ln())
            .collect::<Vec<_>>()
            .into_boxed_slice()
    })
}

/// `log(e^a + e^b)` (natural log) via the table `tbl` (from [`table`]). `-inf` is the additive
/// identity. Bit-for-bit independent of how the caller obtained `tbl`.
#[inline]
pub fn flogsum(tbl: &[f32], a: f32, b: f32) -> f32 {
    if a == f32::NEG_INFINITY {
        return b;
    }
    if b == f32::NEG_INFINITY {
        return a;
    }
    let (max, min) = if a >= b { (a, b) } else { (b, a) };
    let d = max - min;
    if d >= MAXDIFF {
        max
    } else {
        // d < MAXDIFF < 16 ⇒ index < 15_700 < TBL, so the lookup is always in bounds.
        max + tbl[(d * SCALE) as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Table lookup matches the exact transcendental within HMMER's documented tolerance.
    #[test]
    fn matches_exact_flogsum() {
        let tbl = table();
        let exact = |a: f32, b: f32| -> f32 {
            let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
            hi + ((lo - hi) as f64).exp().ln_1p() as f32
        };
        // Sweep a range of (a, b) including near-equal and far-apart pairs.
        for ai in -50..50 {
            for bi in -50..50 {
                let a = ai as f32 * 0.37;
                let b = bi as f32 * 0.41;
                let got = flogsum(tbl, a, b);
                let want = exact(a, b);
                assert!(
                    (got - want).abs() < 1e-3,
                    "flogsum({a}, {b}) = {got}, exact {want}"
                );
            }
        }
    }

    #[test]
    fn neg_inf_is_identity() {
        let tbl = table();
        assert_eq!(flogsum(tbl, f32::NEG_INFINITY, 3.5), 3.5);
        assert_eq!(flogsum(tbl, -2.0, f32::NEG_INFINITY), -2.0);
    }
}
