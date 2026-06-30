//! CP9 (CM Plan 9) HMM — the profile-HMM Infernal derives from a CM to drive **HMM banding**
//! (per-sequence, posterior-derived j+d bands that tighten the CYK/Inside DP far beyond the
//! model-only d-bands of [`crate::qdb`]). This is a multi-stage port of Infernal's `cp9*.c` /
//! `hmmband.c`, built up in milestones:
//!
//!   1. **expected state occupancy** (`psi`) — the foundation of the CM→CP9 parameter mapping.  ✓
//!   2. the CP9 HMM construction: (a) the **CM↔CP9 map** (which CM states map to each CP9 node's
//!      match/insert/delete) ← here, then (b) the occupancy-weighted emission + transition
//!      parameters, validated psi-vs-phi.
//!   3. CP9 Forward / Backward + posteriors.
//!   4. posteriors → HMM position bands → CM i,j bands → per-j d-bands (`hdmin`/`hdmax`).
//!   5. the doubly-ragged (j- and d-banded) CYK kernel that consumes them.
//!
//! The CP9 build (1–2) depends only on the CM, so it runs once per model; 3–4 are per sequence.
//!
//! Milestone 1 (this file) is foundation scaffolding for the later stages; its items are
//! exercised by tests now and consumed by the CP9 construction next, hence `allow(dead_code)`.
#![allow(dead_code)]

use crate::emitmap::EmitMap;
use crate::model::{nd, st, Cm};

/// CP9 node-internal state types (index into `hns2cs[k][_]` and the CP9 transition tables).
pub(crate) const HMM_MATCH: usize = 0;
pub(crate) const HMM_INSERT: usize = 1;
pub(crate) const HMM_DELETE: usize = 2;

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

/// The CM↔CP9 correspondence (`CP9Map_t`). A CP9 HMM has one node per CM **consensus column**
/// (`clen` nodes, plus node 0 for the ROOT inserts); each node has a match / insert / delete
/// state. A CM MATL/MATR node contributes one column, a MATP node two (left = `lpos`, right =
/// `rpos`). The maps are bidirectional so the parameter build can pull CM occupancy/emissions
/// into CP9 states and validate the reverse.
pub(crate) struct Cp9Map {
    /// Number of CP9 nodes = consensus length (`clen`).
    pub hmm_m: usize,
    /// `[cm_node]` left / right consensus column for MAT* subtrees, else `-1`.
    pub nd2lpos: Vec<i32>,
    pub nd2rpos: Vec<i32>,
    /// `[1..=clen]` the MATL/MATR/MATP CM node modelling this consensus column.
    pub pos2nd: Vec<i32>,
    /// `[cm_state][0..2]` the 1–2 CP9 nodes a CM state maps to (`[_][1] == -1` if just one).
    pub cs2hn: Vec<[i32; 2]>,
    /// `[cm_state][0..2]` the matching CP9 state types (`HMM_MATCH|INSERT|DELETE`).
    pub cs2hs: Vec<[i32; 2]>,
    /// `[hmm_node][state_type][0..2]` the 1–2 CM states mapping to that CP9 node-state.
    pub hns2cs: Vec<[[i32; 2]; 3]>,
}

/// Build the CM↔CP9 map — a port of `CP9_map_cm2hmm`. Consensus columns come from the emit map;
/// the per-node-type casework records, for each CP9 node, which CM split-set state is its match,
/// which CM insert reaches "between" the columns, and which CM deletes/half-matches it covers.
///
/// The fiddly part is insert attachment: a CM MATP_IR / MATR_IR inserts to the *left* of its
/// column, but a CP9 node's insert emits to the *right*, so an IR is attached to the CP9 node of
/// the neighbouring column (found by inspecting which CM node models column `k+1`).
pub(crate) fn build_cp9_map(cm: &Cm) -> Cp9Map {
    let clen = cm.clen;
    let emap = EmitMap::build(cm);
    let mut nd2lpos = vec![-1i32; cm.nodes];
    let mut nd2rpos = vec![-1i32; cm.nodes];
    let mut pos2nd = vec![-1i32; clen + 1];
    for n in 0..cm.nodes {
        if cm.ndtype[n] == nd::MATP || cm.ndtype[n] == nd::MATL {
            nd2lpos[n] = emap.lpos[n];
            pos2nd[emap.lpos[n] as usize] = n as i32;
        }
        if cm.ndtype[n] == nd::MATP || cm.ndtype[n] == nd::MATR {
            nd2rpos[n] = emap.rpos[n];
            pos2nd[emap.rpos[n] as usize] = n as i32;
        }
    }

    let mut map = Cp9Map {
        hmm_m: clen,
        nd2lpos,
        nd2rpos,
        pos2nd,
        cs2hn: vec![[-1, -1]; cm.m],
        cs2hs: vec![[-1, -1]; cm.m],
        hns2cs: vec![[[-1, -1]; 3]; clen + 1],
    };

    // Records the CM state v as CP9 node k's state-type `ks` (and the reverse). An insert state
    // immediately before an END is "detached" (removed to break score ambiguities) and skipped.
    let mut map_helper = |cm: &Cm, m: &mut Cp9Map, k: usize, ks: usize, v: usize| {
        if ks == HMM_INSERT && cm.sttype[v + 1] == st::E {
            return;
        }
        if m.cs2hn[v][0] == -1 {
            m.cs2hn[v][0] = k as i32;
            m.cs2hs[v][0] = ks as i32;
        } else if m.cs2hn[v][1] == -1 {
            m.cs2hn[v][1] = k as i32;
            m.cs2hs[v][1] = ks as i32;
        } else {
            panic!("cp9map: cs2hn[{v}] overflow");
        }
        if m.hns2cs[k][ks][0] == -1 {
            m.hns2cs[k][ks][0] = v as i32;
        } else if m.hns2cs[k][ks][1] == -1 {
            m.hns2cs[k][ks][1] = v as i32;
        } else {
            panic!("cp9map: hns2cs[{k}][{ks}] overflow");
        }
    };

    // Node 0: ROOT_S → match, ROOT_IL → insert. (ROOT_IR is attached to the last node below.)
    map_helper(cm, &mut map, 0, HMM_MATCH, 0);
    map_helper(cm, &mut map, 0, HMM_INSERT, 1);

    let nm = |n: usize| cm.nodemap[n]; // first CM state of node n
    for k in 1..=map.hmm_m {
        let n = map.pos2nd[k] as usize;
        let is_left = map.nd2lpos[n] == k as i32;
        match cm.ndtype[n] {
            nd::MATP if is_left => {
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n)); // MATP_MP
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n) + 1); // MATP_ML
                map_helper(cm, &mut map, k, HMM_INSERT, nm(n) + 4); // MATP_IL
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 2); // MATP_MR
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 3); // MATP_D
            }
            nd::MATP => {
                // right half of a MATP
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n)); // MATP_MP
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n) + 2); // MATP_MR
                attach_right_insert(cm, &mut map, k, &mut map_helper);
                // The MATP_IR also belongs to the previous column's CP9 insert if that column is
                // a non-MATL left half (otherwise its contribution would be lost).
                let pk1 = map.pos2nd[k - 1] as usize;
                if map.nd2lpos[pk1] == (k as i32 - 1) && cm.ndtype[pk1] != nd::MATL {
                    assert_eq!(cm.ndtype[pk1], nd::MATP, "cp9map: unexpected k-1 node");
                    map_helper(cm, &mut map, k - 1, HMM_INSERT, nm(n) + 5); // MATP_IR
                }
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 1); // MATP_ML
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 3); // MATP_D
            }
            nd::MATL => {
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n)); // MATL_ML
                map_helper(cm, &mut map, k, HMM_INSERT, nm(n) + 2); // MATL_IL
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 1); // MATL_D
                if k == map.hmm_m {
                    map_helper(cm, &mut map, k, HMM_INSERT, 2); // ROOT_IR
                }
            }
            nd::MATR => {
                map_helper(cm, &mut map, k, HMM_MATCH, nm(n)); // MATR_MR
                attach_right_insert(cm, &mut map, k, &mut map_helper);
                map_helper(cm, &mut map, k, HMM_DELETE, nm(n) + 1); // MATR_D
            }
            other => panic!("cp9map: HMM node {k} maps to non-MAT node type {other}"),
        }
    }
    map
}

/// Attach the CP9 *insert* of node `k` (which emits to the right of column `k`). The CM state
/// that inserts there is found from the node modelling column `k+1`: if `k+1` is a left half,
/// it's the nearest enclosing BEGR's IL; if a right half, it's that node's IR. At the last node
/// the right insert is ROOT_IR. Shared by the MATP-right and MATR cases.
fn attach_right_insert(
    cm: &Cm,
    map: &mut Cp9Map,
    k: usize,
    map_helper: &mut impl FnMut(&Cm, &mut Cp9Map, usize, usize, usize),
) {
    if k == map.hmm_m {
        map_helper(cm, map, k, HMM_INSERT, 2); // ROOT_IR
        return;
    }
    let nn = map.pos2nd[k + 1] as usize;
    if map.nd2lpos[nn] == (k as i32 + 1) {
        // column k+1 is a left half: the insert is the IL of the nearest BEGR at/above nn.
        let mut nbegr = nn as i32;
        while nbegr >= 0 && cm.ndtype[nbegr as usize] != nd::BEGR {
            nbegr -= 1;
        }
        assert!(nbegr >= 0, "cp9map: no BEGR above node {nn}");
        map_helper(cm, map, k, HMM_INSERT, cm.nodemap[nbegr as usize] + 1); // BEGR_IL
    } else if map.nd2rpos[nn] == (k as i32 + 1) {
        match cm.ndtype[nn] {
            nd::MATP => map_helper(cm, map, k, HMM_INSERT, cm.nodemap[nn] + 5), // MATP_IR
            nd::MATR => map_helper(cm, map, k, HMM_INSERT, cm.nodemap[nn] + 2), // MATR_IR
            _ => {}
        }
    }
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

    /// The CM↔CP9 map must be complete and self-consistent: every consensus column maps to a
    /// MAT* node, every CP9 node has a match state mapping to a CM consensus emitter, and the
    /// forward (cs2h*) and reverse (hns2cs) maps agree.
    #[test]
    fn cp9_map_is_consistent() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let map = build_cp9_map(&cm);

        assert_eq!(map.hmm_m, cm.clen, "one CP9 node per consensus column");

        for k in 1..=map.hmm_m {
            // Each consensus column is modelled by a MATL/MATR/MATP node.
            let n = map.pos2nd[k];
            assert!(n >= 0, "column {k} maps to no CM node");
            let ndt = cm.ndtype[n as usize];
            assert!(
                ndt == nd::MATP || ndt == nd::MATL || ndt == nd::MATR,
                "column {k} maps to non-MAT node type {ndt}"
            );
            // Each CP9 node has a match state, and it maps to a CM consensus emitter (MP/ML/MR).
            let m0 = map.hns2cs[k][HMM_MATCH][0];
            assert!(m0 >= 0, "CP9 node {k} has no match state");
            let mt = cm.sttype[m0 as usize];
            assert!(
                mt == st::MP || mt == st::ML || mt == st::MR,
                "CP9 node {k} match maps to non-emitter state type {mt}"
            );
        }

        // Bidirectional consistency: every CM→CP9 edge appears in the reverse map.
        for v in 0..cm.m {
            for slot in 0..2 {
                let k = map.cs2hn[v][slot];
                if k < 0 {
                    continue;
                }
                let ks = map.cs2hs[v][slot];
                assert!(
                    ks >= 0,
                    "state {v} has node but no state-type in slot {slot}"
                );
                let (k, ks) = (k as usize, ks as usize);
                let rev = map.hns2cs[k][ks];
                assert!(
                    rev[0] == v as i32 || rev[1] == v as i32,
                    "reverse map hns2cs[{k}][{ks}] missing CM state {v}"
                );
            }
        }
    }
}
