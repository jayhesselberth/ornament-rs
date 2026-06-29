//! `cm_Configure` equivalent: derive log-odds score vectors from the parsed
//! probabilities. Faithful port of the score computation in `CMLogoddsify` (`src/cm.c`).
//!
//! The `.cm` file stores log-odds scores; the parser converts them to probabilities, and
//! this step converts them back to scores — the round trip reproduces the file values
//! exactly (to float precision). Local/glocal begin/end configuration and the degenerate
//! emission expansion (`oesc`) are layered on in later phases; this provides the canonical
//! `tsc`/`esc` the scanning DP consumes.

use crate::model::{nd, st, Cm};

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

/// Configure the model for **local** alignment (`cm_localize` + `cm_CalculateLocalBeginProbs`).
///
/// Spreads `pbegin` local-entry mass across internal MATP/MATL/MATR/BIF nodes and `pend`
/// local-exit mass across internal MATP/MATL/MATR/BEGL/BEGR nodes (renormalizing their
/// transitions), zeroes the ROOT transitions (entry is only via local begins), then
/// recomputes all log-odds scores. Call this INSTEAD of [`configure_scores`] for local mode.
pub fn configure_local(cm: &mut Cm) {
    if cm.is_local {
        return;
    }
    let p_start = cm.pbegin;
    let p_exit = cm.pend;
    let is_start_nd = |t: u8| matches!(t, nd::MATP | nd::MATL | nd::MATR | nd::BIF);
    let is_exit_nd = |t: u8| matches!(t, nd::MATP | nd::MATL | nd::MATR | nd::BEGL | nd::BEGR);

    // --- Local begins: spread p_start over internal start nodes (nd >= 2); node 1 gets the rest.
    let nstarts = (2..cm.nodes).filter(|&nd_i| is_start_nd(cm.ndtype[nd_i])).count();
    for b in cm.begin.iter_mut() {
        *b = 0.0;
    }
    if cm.nodes > 1 {
        cm.begin[cm.nodemap[1]] = 1.0 - p_start;
    }
    if nstarts > 0 {
        let p = p_start / nstarts as f32;
        for nd_i in 2..cm.nodes {
            if is_start_nd(cm.ndtype[nd_i]) {
                cm.begin[cm.nodemap[nd_i]] = p;
            }
        }
    }
    // Entry to ROOT's children is removed; the only way out of ROOT is a local begin.
    for x in cm.t[0].iter_mut() {
        *x = 0.0;
    }

    // --- Local ends: spread p_exit over internal exit nodes not adjacent to END.
    let nexits = (1..cm.nodes)
        .filter(|&nd_i| {
            is_exit_nd(cm.ndtype[nd_i]) && nd_i + 1 < cm.nodes && cm.ndtype[nd_i + 1] != nd::END
        })
        .count();
    for e in cm.end.iter_mut() {
        *e = 0.0;
    }
    if nexits > 0 {
        for nd_i in 1..cm.nodes {
            if is_exit_nd(cm.ndtype[nd_i]) && nd_i + 1 < cm.nodes && cm.ndtype[nd_i + 1] != nd::END {
                let v = cm.nodemap[nd_i];
                cm.end[v] = p_exit / nexits as f32;
                let denom: f32 = cm.t[v].iter().sum::<f32>() + cm.end[v];
                if denom > 0.0 {
                    for x in cm.t[v].iter_mut() {
                        *x /= denom;
                    }
                }
            }
        }
    }

    // Ensure el_selfsc * W can't underflow past IMPOSSIBLE (matches cm_Configure).
    if (cm.el_selfsc * cm.w as f32) < IMPOSSIBLE {
        cm.el_selfsc = IMPOSSIBLE / (cm.w as f32 + 1.0);
    }

    cm.is_local = true;

    // Recompute tsc/esc from the modified probabilities, plus the begin/end scores.
    configure_scores(cm);
    for v in 0..cm.m {
        cm.beginsc[v] = sre_log2(cm.begin[v]);
        cm.endsc[v] = sre_log2(cm.end[v]);
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
