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
//! These items are foundation scaffolding for the later stages — exercised by tests now and
//! consumed by the CP9 construction next, hence `allow(dead_code)`. Residue-indexed loops over
//! the alphabet read clearer as range loops (the index feeds `emit_prob` and indexes the array).
#![allow(dead_code, clippy::needless_range_loop)]

use crate::emitmap::EmitMap;
use crate::model::{nd, st, Cm};

/// CP9 node-internal state types (index into `hns2cs[k][_]` and the CP9 transition tables).
pub(crate) const HMM_MATCH: usize = 0;
pub(crate) const HMM_INSERT: usize = 1;
pub(crate) const HMM_DELETE: usize = 2;

/// CP9 transition indices into `Cp9::t[k]` (`CT**`). Match block = `[0..4]` (+ `end[k]`),
/// insert block = `[4..7]`, delete block = `[7..10]`.
pub(crate) const CTMM: usize = 0; // M_k -> M_{k+1}
pub(crate) const CTMI: usize = 1; // M_k -> I_k
pub(crate) const CTMD: usize = 2; // M_k -> D_{k+1}
pub(crate) const CTMEL: usize = 3; // M_k -> EL (local end; 0 in glocal)
pub(crate) const CTIM: usize = 4; // I_k -> M_{k+1}
pub(crate) const CTII: usize = 5; // I_k -> I_k
pub(crate) const CTID: usize = 6; // I_k -> D_{k+1}
pub(crate) const CTDM: usize = 7; // D_k -> M_{k+1}
pub(crate) const CTDI: usize = 8; // D_k -> I_k
pub(crate) const CTDD: usize = 9; // D_k -> D_{k+1}

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
                // The MATP_IR inserts between columns k-1 and k. Whenever k-1 is a left-half
                // column (MATL or MATP-left), the IR attaches to that previous CP9 node's insert
                // — including the MATL case, where the column's own (detached) MATL_IL can't model
                // this insertion, so the MATP_IR is what fills the node's insert.
                let pk1 = map.pos2nd[k - 1] as usize;
                if map.nd2lpos[pk1] == (k as i32 - 1) {
                    debug_assert!(cm.ndtype[pk1] == nd::MATL || cm.ndtype[pk1] == nd::MATP);
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

/// A CM Plan 9 HMM (`CP9_t`), in probability form. `m` nodes (one per consensus column); node 0
/// carries only the leading insert. Match/insert emissions are over the `k`-symbol alphabet; the
/// nine Plan-9 transitions per node and the begin/end distributions are filled by the transition
/// build (milestone 2b-ii) and are left at zero here.
pub(crate) struct Cp9 {
    pub m: usize,
    pub k: usize,
    /// `[0..=m][0..k]` match emissions (`mat[0]` unused).
    pub mat: Vec<Vec<f64>>,
    /// `[0..=m][0..k]` insert emissions (`ins[0]` = the leading ROOT_IL insert).
    pub ins: Vec<Vec<f64>>,
    /// `[0..=m][0..10]` Plan-9 transitions, indexed by the `CT**` constants.
    pub t: Vec<[f64; 10]>,
    /// `[1..=m]` B→Mk begin / Mk→E end distributions (begin/last-node end in 2b-ii special).
    pub begin: Vec<f64>,
    pub end: Vec<f64>,
}

/// Emission probability of residue `i` from CM state `x` as seen by CP9 node `k` — a port of
/// `cm2hmm_emit_prob`. Singlet emitters pass through `cm.e`; a pair emitter (MATP_MP) is
/// **marginalized** onto the column the node models: the left column sums over the right base
/// (`e[x][i*K .. (i+1)*K]`), the right column strides over the left base (`e[x][i], e[x][i+K], …`).
fn emit_prob(cm: &Cm, map: &Cp9Map, x: usize, i: usize, k: usize) -> f64 {
    if cm.stid[x] != crate::model::stid::MATP_MP {
        return cm.e[x][i] as f64;
    }
    let kk = cm.k();
    let is_left = map.nd2lpos[map.pos2nd[k] as usize] == k as i32;
    if is_left {
        (i * kk..(i + 1) * kk).map(|j| cm.e[x][j] as f64).sum()
    } else {
        (i..kk * kk).step_by(kk).map(|j| cm.e[x][j] as f64).sum()
    }
}

/// Build the CP9 **emissions** (milestone 2b-i): each CP9 match/insert emission is the
/// occupancy-weighted average of the CM emissions that map to it (`Σ psi[ap]·emit_prob(ap,·)`),
/// then renormalized to a distribution (the `CPlan9Renormalize` step for emissions). The leading
/// insert (node 0) is special-cased to ROOT_IL (CM state 1), per `build_cp9_hmm`. Transitions are
/// left at zero for the next milestone.
pub(crate) fn build_cp9_emissions(cm: &Cm) -> Cp9 {
    let map = build_cp9_map(cm);
    let psi = expected_occupancy(cm);
    let k = cm.k();
    let m = map.hmm_m;

    let mut mat = vec![vec![0.0f64; k]; m + 1];
    let mut ins = vec![vec![0.0f64; k]; m + 1];

    // Node 0's insert maps to ROOT_IL = CM state 1.
    for i in 0..k {
        ins[0][i] = cm.e[1][i] as f64;
    }
    for node in 1..=m {
        for slot in 0..2 {
            let ap = map.hns2cs[node][HMM_MATCH][slot];
            if ap >= 0 {
                for i in 0..k {
                    mat[node][i] += psi[ap as usize] * emit_prob(cm, &map, ap as usize, i, node);
                }
            }
            let ai = map.hns2cs[node][HMM_INSERT][slot];
            if ai >= 0 {
                for i in 0..k {
                    ins[node][i] += psi[ai as usize] * emit_prob(cm, &map, ai as usize, i, node);
                }
            }
        }
    }

    // Renormalize each emission vector to a probability distribution.
    let renorm = |v: &mut [f64]| {
        let s: f64 = v.iter().sum();
        if s > 0.0 {
            v.iter_mut().for_each(|x| *x /= s);
        }
    };
    for node in 0..=m {
        renorm(&mut ins[node]);
        if node >= 1 {
            renorm(&mut mat[node]);
        }
    }

    Cp9 {
        m,
        k,
        mat,
        ins,
        t: vec![[0.0; 10]; m + 1],
        begin: vec![0.0; m + 1],
        end: vec![0.0; m + 1],
    }
}

/// Summed probability of all CM paths from `start` to `end` (`start ≤ end`, i.e. *down* the CM)
/// — a port of `cm_sum_subpaths_cp9`. It runs a localized occupancy recursion over the states in
/// `[start, end]`: `sub[0]=1` at `start`, each later state accumulates `Σ sub[parent]·t[parent→·]`
/// (parents before `start` are out of window and contribute 0), inserts fold their self-loop.
///
/// Two corrections keep this consistent with the way the 9 CP9 transitions partition CM paths:
///  - when both ends are non-inserts, the mass that reaches `start` *via this node's own insert*
///    is subtracted (it belongs to the M→I / D→I transitions, counted separately);
///  - an interior insert state that maps to node `k` is skipped (counted in its own sub-call).
///
/// The two degenerate `start == end` cases: a self-insert returns its self-loop prob; a non-insert
/// (only happens for the two halves of one MATP modelling adjacent columns) returns 1.
fn sum_subpaths(cm: &Cm, map: &Cp9Map, start: usize, end: usize, k: usize, psi: &[f64]) -> f64 {
    debug_assert!(start <= end);
    let is_ins = |v: usize| cm.sttype[v] == st::IL || cm.sttype[v] == st::IR;

    if start == end {
        if !is_ins(start) {
            return 1.0; // adjacent columns share a MATP node; caller weights by psi[start]
        }
        return cm.t[start][0] as f64; // self-insert probability
    }

    let mut sub = vec![0.0f64; end - start + 1];
    sub[0] = 1.0; // we begin in `start`

    // Subtract paths that arrive at `start` through this node's insert (counted by M→I / D→I).
    if !is_ins(start) && !is_ins(end) {
        let mut insert_to_start = 0.0;
        for slot in 0..2 {
            let iv = map.hns2cs[k][HMM_INSERT][slot];
            if iv >= 0 && (iv as usize) < start {
                let iv = iv as usize;
                insert_to_start += psi[iv] * sum_subpaths(cm, map, iv, start, k, psi);
            }
        }
        sub[0] -= insert_to_start / psi[start];
    }

    for v in (start + 1)..=end {
        let isv = is_ins(v);
        if cm.sttype[v] == st::S {
            // A start state is entered structurally from a B (prob 1) — its only "parent" is that
            // bifurcation, which carries no transition probabilities (cm.t[B] is empty), so the
            // parent loop would contribute 0. Carry the prior value, matching `sub_psi[v-1]*1`.
            sub[v - start] = sub[v - 1 - start];
            continue;
        }
        // Skip an interior node-k insert: its contribution is counted in its own sub-call.
        let skip = v != end && isv && (map.cs2hn[v][0] == k as i32 || map.cs2hn[v][1] == k as i32);
        if skip {
            continue;
        }
        let final_y = if isv { 1 } else { 0 };
        for y in final_y..cm.pnum[v] {
            let x = (cm.plast[v] - y) as usize;
            if x >= start {
                let off = v - cm.cfirst[x] as usize;
                sub[v - start] += sub[x - start] * cm.t[x][off] as f64;
            }
        }
        if v != end && isv {
            let ts = cm.t[v][0] as f64; // self-loop, folded in closed form
            sub[v - start] += sub[v - start] * (ts / (1.0 - ts));
        }
    }
    sub[end - start]
}

/// Add the virtual count for one CP9 transition `a → b` (CM states) — a port of
/// `hmm_add_single_trans_cp9`. The path always runs down the CM (low→high index), weighted by the
/// occupancy of its upstream end.
#[allow(clippy::too_many_arguments)]
fn add_single_trans(
    cm: &Cm,
    map: &Cp9Map,
    hmm: &mut Cp9,
    a: i32,
    b: i32,
    k: usize,
    idx: usize,
    psi: &[f64],
) {
    if a < 0 || b < 0 {
        return;
    }
    let (a, b) = (a as usize, b as usize);
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    hmm.t[k][idx] += psi[lo] * sum_subpaths(cm, map, lo, hi, k, psi);
}

/// The 9 interior CP9 transitions of node `k`, as `(a-state-type, b-state-type, b-node-offset)`.
const TRANS_SPEC: [(usize, usize, usize, usize); 9] = [
    (HMM_MATCH, HMM_MATCH, 1, CTMM),
    (HMM_MATCH, HMM_INSERT, 0, CTMI),
    (HMM_MATCH, HMM_DELETE, 1, CTMD),
    (HMM_INSERT, HMM_MATCH, 1, CTIM),
    (HMM_INSERT, HMM_INSERT, 0, CTII),
    (HMM_INSERT, HMM_DELETE, 1, CTID),
    (HMM_DELETE, HMM_MATCH, 1, CTDM),
    (HMM_DELETE, HMM_INSERT, 0, CTDI),
    (HMM_DELETE, HMM_DELETE, 1, CTDD),
];

/// Fill and normalize the 9 transitions of an interior node `k` (`1 ≤ k < M`) — a port of
/// `cm2hmm_trans_probs_cp9`. Each transition's virtual count sums over the (≤2)×(≤2) CM-state
/// pairs its endpoints map to; the three out-blocks (match incl. `end[k]`, insert, delete) are
/// then renormalized to distributions. Node 0 and node M are handled by the special build.
fn fill_interior_transitions(cm: &Cm, map: &Cp9Map, hmm: &mut Cp9, k: usize, psi: &[f64]) {
    for &(a_ty, b_ty, b_off, idx) in TRANS_SPEC.iter() {
        let ap = map.hns2cs[k][a_ty];
        let bp = map.hns2cs[k + b_off][b_ty];
        for &a in ap.iter() {
            for &b in bp.iter() {
                add_single_trans(cm, map, hmm, a, b, k, idx, psi);
            }
        }
    }
    // Renormalize the three out-of-state blocks.
    let d = hmm.t[k][0..4].iter().sum::<f64>() + hmm.end[k];
    if d > 0.0 {
        for x in hmm.t[k][0..4].iter_mut() {
            *x /= d;
        }
        hmm.end[k] /= d;
    }
    for block in [4usize, 7] {
        let s: f64 = hmm.t[k][block..block + 3].iter().sum();
        if s > 0.0 {
            for x in hmm.t[k][block..block + 3].iter_mut() {
                *x /= s;
            }
        }
    }
}

/// Build the CP9 emissions and the **interior-node transitions** (`1 ≤ k < M`). Node 0 / node M
/// special transitions and the begin/end distributions are the next milestone, so the boundary
/// nodes' transitions are left zero here.
pub(crate) fn build_cp9_interior(cm: &Cm) -> Cp9 {
    let map = build_cp9_map(cm);
    let psi = expected_occupancy(cm);
    let mut hmm = build_cp9_emissions(cm);
    for k in 1..hmm.m {
        fill_interior_transitions(cm, &map, &mut hmm, k, &psi);
    }
    hmm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmfile::parse_cm_file;
    use crate::config::configure_scores;
    use crate::model::stid;

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

    /// CP9 emissions must be probability distributions, and a singlet (MATL) column must
    /// reproduce the CM's consensus emission exactly (its match maps to one CM state, MATL_ML).
    #[test]
    fn cp9_emissions_are_distributions() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let map = build_cp9_map(&cm);
        let hmm = build_cp9_emissions(&cm);

        for node in 1..=hmm.m {
            let s: f64 = hmm.mat[node].iter().sum();
            assert!((s - 1.0).abs() < 1e-9, "mat[{node}] sums to {s}");
        }
        for node in 0..=hmm.m {
            let s: f64 = hmm.ins[node].iter().sum();
            assert!((s - 1.0).abs() < 1e-9, "ins[{node}] sums to {s}");
        }

        // A MATL column (match maps to a single MATL_ML state) reproduces the CM emission.
        for node in 0..=hmm.m {
            if hmm.ins[node].iter().sum::<f64>() == 0.0 {
                let n = map.pos2nd[node];
                let nt = if n >= 0 {
                    cm.ndtype[n as usize] as i32
                } else {
                    -1
                };
                let nn = if node < hmm.m {
                    map.pos2nd[node + 1]
                } else {
                    -1
                };
                let nnt = if nn >= 0 {
                    cm.ndtype[nn as usize] as i32
                } else {
                    -1
                };
                eprintln!("DBG empty-ins node={node} ndtype={nt} hns2cs_ins={:?} k+1_nd={nn} k+1_type={nnt} nd2lpos[k+1nd]={} nd2rpos[k+1nd]={}",
                    map.hns2cs[node][HMM_INSERT], if nn>=0 {map.nd2lpos[nn as usize]} else {-99}, if nn>=0 {map.nd2rpos[nn as usize]} else {-99});
            }
        }
        let mut checked = 0;
        for node in 1..=hmm.m {
            let n = map.pos2nd[node] as usize;
            if cm.ndtype[n] == nd::MATL && map.hns2cs[node][HMM_MATCH][1] == -1 {
                let ml = cm.nodemap[n]; // MATL_ML
                assert_eq!(cm.stid[ml], stid::MATL_ML);
                // Compare against the *renormalized* CM emission (the CP9 build renormalizes, and
                // the stored f32 distribution doesn't sum to exactly 1).
                let denom: f64 = (0..hmm.k).map(|i| cm.e[ml][i] as f64).sum();
                for i in 0..hmm.k {
                    let want = cm.e[ml][i] as f64 / denom;
                    assert!(
                        (hmm.mat[node][i] - want).abs() < 1e-9,
                        "MATL column {node} emission[{i}] {} != cm.e {}",
                        hmm.mat[node][i],
                        want
                    );
                }
                checked += 1;
            }
        }
        assert!(checked > 0, "expected at least one MATL column to validate");
    }

    /// Interior-node transitions must form valid distributions: the match out-block (`[0..4]` +
    /// `end[k]`), the insert out-block (`[4..7]`), and the delete out-block (`[7..10]`) each sum
    /// to 1. Also sanity-check that a confident column prefers M→M.
    #[test]
    fn cp9_interior_transitions_normalize() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let hmm = build_cp9_interior(&cm);

        let mut mm_dominant = 0;
        for k in 1..hmm.m {
            let mblock: f64 = hmm.t[k][0..4].iter().sum::<f64>() + hmm.end[k];
            assert!(
                (mblock - 1.0).abs() < 1e-9,
                "node {k} match block = {mblock}"
            );
            let iblock: f64 = hmm.t[k][4..7].iter().sum();
            assert!(
                (iblock - 1.0).abs() < 1e-9,
                "node {k} insert block = {iblock}"
            );
            let dblock: f64 = hmm.t[k][7..10].iter().sum();
            assert!(
                (dblock - 1.0).abs() < 1e-9,
                "node {k} delete block = {dblock}"
            );
            // Probabilities are in [0, 1].
            for &x in hmm.t[k].iter() {
                assert!(
                    (0.0..=1.0).contains(&x),
                    "node {k} transition {x} out of range"
                );
            }
            if hmm.t[k][CTMM] > hmm.t[k][CTMI] && hmm.t[k][CTMM] > hmm.t[k][CTMD] {
                mm_dominant += 1;
            }
        }
        // A well-trained consensus model: most columns most-prefer staying in the match path.
        assert!(
            mm_dominant > hmm.m / 2,
            "expected M->M dominant in most columns, got {mm_dominant}/{}",
            hmm.m - 1
        );
    }
}
