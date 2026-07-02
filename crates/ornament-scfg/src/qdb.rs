//! Query-dependent banding (QDB) — per-state subsequence-length bands `[dmin[v], dmax[v]]`
//! that restrict the CYK/Inside `d`-loops to the plausible length range, turning the scan's
//! inner cost from O(L·W²·M) toward O(L·W·band·M). Port of `BandCalculationEngine`
//! (`cm_qdband.c`).
//!
//! The band for state `v` is derived from `gamma[v][n]` = P(a parse subtree rooted at `v`
//! emits a subsequence of length `n`), computed bottom-up over states:
//! - `E`: emits length 0, so `gamma[v][0] = 1`.
//! - `B` (bifurcation): the convolution of its two children's length densities.
//! - everything else: shift the child density by `StateDelta` (MP=2, M*/I*=1, else 0),
//!   weighted by the (local-ends-off, renormalized) transition probabilities.
//!
//! `dmin[v]` / `dmax[v]` are then the smallest/largest `n` at which each cumulative tail still
//! holds more than `beta` probability mass; tail loss `≤ beta` per side. With a safe `beta`
//! the optimal parse of any real hit lies entirely within the bands, so a banded scan reports
//! exactly the same hits as the unbanded one (verified by the parity tests).

use serde::{Deserialize, Serialize};

use crate::config::IMPOSSIBLE;
use crate::model::{n_emit, st, Cm};

/// Per-state length bands plus the overall max hit width.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QdbBands {
    /// `dmin[v]` = minimum subsequence length allowed for state `v`.
    pub dmin: Vec<usize>,
    /// `dmax[v]` = maximum subsequence length allowed for state `v` (clamped to `w`).
    pub dmax: Vec<usize>,
    /// Maximum hit width (the model's `W`); the band-calc length ceiling.
    pub w: usize,
}

impl QdbBands {
    /// The default QDB tail-loss used by `cmsearch` (`--beta`), `1e-7`.
    pub const DEFAULT_BETA: f64 = 1e-7;
}

/// Compute QDB bands for `cm` at tail-loss `beta`. Bands depend only on the model, so this is
/// done once and reused across every window/sequence. `dmax` is clamped to the model's `W`
/// (the scan never considers longer `d`), which also bounds the band-calc length axis.
#[allow(clippy::needless_range_loop)] // length index `n` addresses several gamma rows / tails
pub fn calc_qdb_bands(cm: &Cm, beta: f64) -> QdbBands {
    let m = cm.m;
    let z = cm.w as usize; // length ceiling = max hit width

    // gamma[v][n] (f64 for tail precision). Built bottom-up; kept for all v so bifurcation
    // children and band derivation can read them.
    let mut gamma: Vec<Vec<f64>> = vec![Vec::new(); m];
    let mut dmin = vec![0usize; m];
    let mut dmax = vec![0usize; m];

    for v in (0..m).rev() {
        if v == 0 {
            // ROOT_S: the whole-hit length range. [0, W] is safe (W is the max hit width) and
            // leaves the cheap ROOT loop unrestricted; internal states carry the speedup.
            dmin[0] = 0;
            dmax[0] = z;
            gamma[0] = vec![0.0; z + 1];
            continue;
        }

        let mut g = vec![0.0f64; z + 1];
        match cm.sttype[v] {
            st::E => g[0] = 1.0,
            st::B => {
                let left = cm.cfirst[v] as usize;
                let right = cm.cnum[v] as usize;
                for n in 0..=z {
                    let mut s = 0.0;
                    for k in 0..=n {
                        s += gamma[left][k] * gamma[right][n - k];
                    }
                    g[n] = s;
                }
            }
            sttype => {
                // Emitter / D / S: emits `sd` residues, then transitions to a child.
                let sd = n_emit(sttype);
                let ych = cm.cfirst[v] as usize;
                let cnum = cm.cnum[v] as usize;
                // Transition probabilities with local ends off: 2^tsc, renormalized to sum 1.
                let mut t = vec![0.0f64; cnum];
                let mut sum = 0.0;
                for (yoff, tp) in t.iter_mut().enumerate() {
                    let s = cm.tsc[v][yoff];
                    *tp = if s <= IMPOSSIBLE {
                        0.0
                    } else {
                        (s as f64).exp2()
                    };
                    sum += *tp;
                }
                if sum > 0.0 {
                    for tp in t.iter_mut() {
                        *tp /= sum;
                    }
                }
                // `n` ascending so a self-transition (insert states: a child == v) reads this
                // same row at the already-finalized `n - sd`, mirroring the DP recurrence.
                for n in sd..=z {
                    let mut s = 0.0;
                    for (yoff, &tp) in t.iter().enumerate() {
                        let child = ych + yoff;
                        let gc = if child == v {
                            g[n - sd]
                        } else {
                            gamma[child][n - sd]
                        };
                        s += tp * gc;
                    }
                    g[n] = s;
                }
            }
        }

        // Left bound: smallest n at which the cumulative low-tail mass first exceeds beta.
        let mut pdf = 0.0;
        for n in 0..=z {
            pdf += g[n];
            if pdf > beta {
                dmin[v] = n;
                break;
            }
        }
        // Right bound: largest n at which the cumulative high-tail mass (summed from the top)
        // still exceeds beta. Summing from the top avoids the 1 - L(n) roundoff.
        let mut pdf = 0.0;
        for n in (0..=z).rev() {
            pdf += g[n];
            if pdf > beta {
                dmax[v] = n;
                break;
            }
        }
        // Degenerate guard: keep the band non-empty and ordered.
        if dmax[v] < dmin[v] {
            dmax[v] = dmin[v];
        }
        gamma[v] = g;
    }

    QdbBands { dmin, dmax, w: z }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmfile::parse_cm_file;
    use crate::config::configure_local;

    #[test]
    fn bands_are_sane() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_local(&mut cm);
        let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);
        let w = cm.w as usize;

        for v in 0..cm.m {
            assert!(
                bands.dmin[v] <= bands.dmax[v],
                "state {v}: dmin {} > dmax {}",
                bands.dmin[v],
                bands.dmax[v]
            );
            assert!(
                bands.dmax[v] <= w,
                "state {v}: dmax {} > W {w}",
                bands.dmax[v]
            );
            if cm.sttype[v] == st::E {
                assert_eq!((bands.dmin[v], bands.dmax[v]), (0, 0), "E state {v}");
            }
        }
        assert_eq!(
            (bands.dmin[0], bands.dmax[0]),
            (0, w),
            "ROOT band is [0, W]"
        );

        // The bands must shrink the work: the summed band width is far below M*W.
        let total: usize = (0..cm.m).map(|v| bands.dmax[v] - bands.dmin[v] + 1).sum();
        assert!(
            total < cm.m * w / 2,
            "bands not tight: total width {total} vs M*W {}",
            cm.m * w
        );
    }
}
