//! CP9 (CM Plan 9) HMM — the profile-HMM Infernal derives from a CM to drive **HMM banding**
//! (per-sequence, posterior-derived j+d bands that tighten the CYK/Inside DP far beyond the
//! model-only d-bands of [`crate::qdb`]). This is a multi-stage port of Infernal's `cp9*.c` /
//! `hmmband.c`, built up in milestones:
//!
//!   1. **expected state occupancy** (`psi`) — the foundation of the CM→CP9 parameter mapping.  ← here
//!   2. the CP9 HMM construction (states / transitions / emissions, validated psi-vs-phi).
//!   3. CP9 Forward / Backward + posteriors.
//!   4. posteriors → HMM position bands → CM i,j bands → per-j d-bands (`hdmin`/`hdmax`).
//!   5. the doubly-ragged (j- and d-banded) CYK kernel that consumes them.
//!
//! The CP9 build (1–2) depends only on the CM, so it runs once per model; 3–4 are per sequence.
//!
//! Milestone 1 (this file) is foundation scaffolding for the later stages; its items are
//! exercised by tests now and consumed by the CP9 construction next, hence `allow(dead_code)`.
#![allow(dead_code)]

use crate::model::{nd, st, Cm};

/// Number of states in a node of the given type (`TotalStatesInNode`). The node's states are
/// `nodemap[nd] .. nodemap[nd] + states_in_node(ndtype[nd])` — which is *not* the same as
/// `nodemap[nd+1]`, because a bifurcation interleaves the BEGL/BEGR subtrees between a node and
/// its tree-order successor.
pub(crate) fn states_in_node(ndtype: u8) -> usize {
    match ndtype {
        nd::BIF => 1,  // B
        nd::MATP => 6, // MP ML MR D IL IR
        nd::MATL => 3, // ML D IL
        nd::MATR => 3, // MR D IR
        nd::BEGL => 1, // S
        nd::BEGR => 2, // S IL
        nd::ROOT => 3, // S IL IR
        nd::END => 1,  // E
        _ => 0,
    }
}

/// Expected state occupancy `psi[v]` = the expected number of times state `v` is entered in a
/// parse — a port of `cm_ExpectedStateOccupancy` (`cm.c`). It propagates occupancy down the
/// guide tree: a start state is entered exactly once (`psi = 1`); every other state accumulates
/// `psi[parent] * P(parent → v)` over its parents, and insert states fold in their self-loop in
/// closed form (`× 1/(1 − t_self)`).
///
/// Transitions are taken with **local ends off** (each state's out-transitions renormalized to
/// sum to 1, removing any EL mass), matching how Infernal computes the bands' occupancy. The
/// invariant that validates it: the `psi` of a node's split-set (non-insert) states sums to 1.
pub(crate) fn expected_occupancy(cm: &Cm) -> Vec<f64> {
    let m = cm.m;

    // Local ends off: per-state transition probabilities renormalized to sum to 1.
    let t: Vec<Vec<f64>> = (0..m)
        .map(|v| {
            let row: Vec<f64> = cm.t[v].iter().map(|&p| p as f64).collect();
            let s: f64 = row.iter().sum();
            if s > 0.0 {
                row.iter().map(|p| p / s).collect()
            } else {
                row
            }
        })
        .collect();

    let mut psi = vec![0.0f64; m];
    for v in 0..m {
        if cm.sttype[v] == st::S {
            psi[v] = 1.0; // start states are entered once in every parse
            continue;
        }
        let is_insert = cm.sttype[v] == st::IL || cm.sttype[v] == st::IR;
        // Parents of v are the contiguous block ending at `plast[v]`: `plast[v] - y` for
        // `y in final_y..pnum[v]`. For inserts, `y = 0` (the self-loop) is excluded and folded
        // in below as a geometric series.
        let final_y = if is_insert { 1 } else { 0 };
        let plast = cm.plast[v];
        for y in final_y..cm.pnum[v] {
            let x = (plast - y) as usize;
            // Transition x → v is `t[x][v - cfirst[x]]` (v's index in x's contiguous children).
            let off = v - cm.cfirst[x] as usize;
            psi[v] += psi[x] * t[x][off];
        }
        if is_insert {
            let tself = t[v][0]; // an insert's first transition is its self-loop
            psi[v] += psi[v] * (tself / (1.0 - tself));
        }
    }
    psi
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmfile::parse_cm_file;
    use crate::config::configure_scores;

    /// Infernal's own correctness check for `psi`: over every node, the occupancy of the
    /// split-set (non-insert) states must sum to 1 — each node is visited once per parse, its
    /// probability mass distributed across that node's match/delete states.
    #[test]
    fn psi_split_set_sums_to_one() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm); // glocal: raw transition probabilities
        let psi = expected_occupancy(&cm);

        for ndi in 0..cm.nodes {
            let lo = cm.nodemap[ndi];
            let hi = lo + states_in_node(cm.ndtype[ndi]);
            let sum: f64 = (lo..hi)
                .filter(|&v| cm.sttype[v] != st::IL && cm.sttype[v] != st::IR)
                .map(|v| psi[v])
                .sum();
            assert!(
                (sum - 1.0).abs() < 1e-3,
                "node {ndi} (type {}) split-set psi = {sum}, expected 1.0",
                cm.ndtype[ndi]
            );
        }
        assert!((psi[0] - 1.0).abs() < 1e-9, "ROOT_S psi = {}", psi[0]);
    }
}
