//! `cm_Configure` equivalent: derive log-odds score vectors from the parsed
//! probabilities. Faithful port of the score computation in `CMLogoddsify` (`src/cm.c`).
//!
//! The `.cm` file stores log-odds scores; the parser converts them to probabilities, and
//! this step converts them back to scores — the round trip reproduces the file values
//! exactly (to float precision). Local/glocal begin/end configuration and the degenerate
//! emission expansion (`oesc`) are layered on in later phases; this provides the canonical
//! `tsc`/`esc` the scanning DP consumes.

use crate::model::{st, Cm};

/// Infernal's "impossible" score sentinel (`IMPOSSIBLE`).
pub const IMPOSSIBLE: f32 = -1e36;

/// `sreLOG2(x) = log2(x)` for `x > 0`, else [`IMPOSSIBLE`]. Uses Infernal's exact constant
/// (`log(x) * 1.44269504`) so scores reproduce the `.cm` file bit-for-bit.
#[inline]
pub fn sre_log2(x: f32) -> f32 {
    if x > 0.0 {
        ((x as f64).ln() * 1.44269504) as f32
    } else {
        IMPOSSIBLE
    }
}

/// Fill `cm.tsc` and `cm.esc` (log-odds scores) from `cm.t`, `cm.e`, and `cm.null`.
pub fn configure_scores(cm: &mut Cm) {
    let k = cm.k();
    for v in 0..cm.m {
        // Transition scores (all states except bifurcation and end).
        if cm.sttype[v] != st::B && cm.sttype[v] != st::E {
            cm.tsc[v] = cm.t[v].iter().map(|&t| sre_log2(t)).collect();
        }

        // Emission scores.
        match cm.sttype[v] {
            st::MP => {
                let mut esc = vec![0.0f32; k * k];
                for x in 0..k {
                    for y in 0..k {
                        let idx = x * k + y;
                        esc[idx] = sre_log2(cm.e[v][idx] / (cm.null[x] * cm.null[y]));
                    }
                }
                cm.esc[v] = esc;
            }
            st::ML | st::MR | st::IL | st::IR => {
                cm.esc[v] = (0..k).map(|x| sre_log2(cm.e[v][x] / cm.null[x])).collect();
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sre_log2_matches_definition() {
        assert!((sre_log2(0.5) - (-1.0)).abs() < 1e-5); // log2(0.5) = -1
        assert!((sre_log2(1.0)).abs() < 1e-6); // log2(1) = 0
        assert_eq!(sre_log2(0.0), IMPOSSIBLE);
    }
}
