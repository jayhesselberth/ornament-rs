//! CP9 (CM Plan 9) HMM — the profile-HMM Infernal derives from a CM to drive **HMM banding**
//! (per-sequence, posterior-derived j+d bands that tighten the CYK/Inside DP far beyond the
//! model-only d-bands of [`crate::qdb`]). This is a multi-stage port of Infernal's `cp9*.c` /
//! `hmmband.c`, built up in milestones:
//!
//!   1. **expected state occupancy** (`psi`) — the foundation of the CM→CP9 parameter mapping.  ✓
//!   2. the CP9 HMM construction — CM↔CP9 map, occupancy-weighted emissions, subpath-summation
//!      transitions + begin/end — validated by the psi-vs-phi occupancy check (`build_cp9`).  ✓
//!   3. CP9 Forward / Backward + posteriors (per sequence), validated F==B + Σ posteriors==1.  ✓
//!   4. posteriors → HMM position bands → CM i,j bands → per-j d-bands (`hdmin`/`hdmax`).  ✓
//!   5. the doubly-ragged (j- and d-banded) CYK kernel that consumes them.  ✓ (alignment kernel;
//!      scanning + multi-hit pipeline wiring is the remaining integration step)
//!
//! The CP9 build (1–2) depends only on the CM, so it runs once per model; 3–4 are per sequence.
//!
//! These items are foundation scaffolding for the later stages — exercised by tests now and
//! consumed by the CP9 construction next, hence `allow(dead_code)`. Residue-indexed loops over
//! the alphabet read clearer as range loops (the index feeds `emit_prob` and indexes the array).
#![allow(dead_code, clippy::needless_range_loop)]

use crate::config::IMPOSSIBLE;
use crate::emit::Oesc;
use crate::emitmap::EmitMap;
use crate::model::{emits_left, emits_right, n_emit, nd, st, stid, Cm};
use ornament_alphabet::Dsq;

/// `ln(p)`, mapping `p ≤ 0` to −∞ (the log-space sentinel).
#[inline]
fn ln(p: f64) -> f64 {
    if p > 0.0 {
        p.ln()
    } else {
        f64::NEG_INFINITY
    }
}

/// Numerically-stable `log(e^a + e^b)`.
#[inline]
fn logsum(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

#[inline]
fn logsum3(a: f64, b: f64, c: f64) -> f64 {
    logsum(logsum(a, b), c)
}

#[inline]
fn logsum4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    logsum(logsum(a, b), logsum(c, d))
}

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

/// Normalize a 3-transition out-block (`insert` at offset 4, `delete` at offset 7) to sum to 1.
fn normalize_block(t: &mut [f64; 10], off: usize) {
    let s: f64 = t[off..off + 3].iter().sum();
    if s > 0.0 {
        for x in t[off..off + 3].iter_mut() {
            *x /= s;
        }
    }
}

/// Fill the **special** boundary transitions — a port of `cm2hmm_special_trans_cp9`. Node 0 is
/// the begin node (B = ROOT_S, N = ROOT_IL): its M→M virtual count becomes `begin[1]` and the
/// rest are the leading N transitions. Node M is the end node: its M/I/D → E counts become the
/// `end[M]` / insert / delete out-transitions. There is no D_0, and several node-M transitions
/// (M→D, I→D, D→D) are illegal and left at 0. Each boundary node's out-blocks are then normalized.
fn fill_special_transitions(cm: &Cm, map: &Cp9Map, hmm: &mut Cp9, psi: &[f64]) {
    let m = hmm.m;

    // ---- Node 0: transitions into node 1 (B and N are ROOT_S / ROOT_IL). ----
    for &(a_ty, b_ty, b_off, idx) in TRANS_SPEC[0..6].iter() {
        let ap = map.hns2cs[0][a_ty];
        let bp = map.hns2cs[b_off][b_ty];
        for &a in ap.iter() {
            for &b in bp.iter() {
                add_single_trans(cm, map, hmm, a, b, 0, idx, psi);
            }
        }
    }
    // B→M1 is the begin distribution, not a t[0][CTMM].
    hmm.begin[1] = hmm.t[0][CTMM];
    hmm.t[0][CTMM] = 0.0;
    // Normalize node 0: all begins + B→N + B→D1 sum to 1; then the N (insert) out-block.
    let d = hmm.begin[1..=m].iter().sum::<f64>() + hmm.t[0][CTMI] + hmm.t[0][CTMD];
    if d > 0.0 {
        for x in hmm.begin[1..=m].iter_mut() {
            *x /= d;
        }
        hmm.t[0][CTMI] /= d;
        hmm.t[0][CTMD] /= d;
    }
    normalize_block(&mut hmm.t[0], 4);

    // ---- Node M: transitions out to the end (E = last CM state). ----
    let e = (cm.m - 1) as i32; // END_E
                               // (a-state, b: None = →E, Some = node-M state, idx)
    let m_pairs = [
        (HMM_MATCH, None, CTMM),
        (HMM_MATCH, Some(HMM_INSERT), CTMI),
        (HMM_INSERT, None, CTIM),
        (HMM_INSERT, Some(HMM_INSERT), CTII),
        (HMM_DELETE, None, CTDM),
        (HMM_DELETE, Some(HMM_INSERT), CTDI),
    ];
    for (a_ty, b_opt, idx) in m_pairs {
        let ap = map.hns2cs[m][a_ty];
        for &a in ap.iter() {
            match b_opt {
                None => add_single_trans(cm, map, hmm, a, e, m, idx, psi),
                Some(b_ty) => {
                    for &b in map.hns2cs[m][b_ty].iter() {
                        add_single_trans(cm, map, hmm, a, b, m, idx, psi);
                    }
                }
            }
        }
    }
    // M→E is the end probability, not t[M][CTMM].
    hmm.end[m] = hmm.t[m][CTMM];
    hmm.t[m][CTMM] = 0.0;
    // Normalize node M's three out-blocks (match incl. end[M], insert, delete).
    let d = hmm.t[m][0..4].iter().sum::<f64>() + hmm.end[m];
    if d > 0.0 {
        for x in hmm.t[m][0..4].iter_mut() {
            *x /= d;
        }
        hmm.end[m] /= d;
    }
    normalize_block(&mut hmm.t[m], 4);
    normalize_block(&mut hmm.t[m], 7);
}

/// Build the complete CP9 HMM from a CM (emissions + all transitions + begin/end) — the
/// per-model parameterization (`build_cp9_hmm`). Runs once per model, not per sequence.
pub(crate) fn build_cp9(cm: &Cm) -> Cp9 {
    let map = build_cp9_map(cm);
    let psi = expected_occupancy(cm);
    let mut hmm = build_cp9_emissions(cm);
    for k in 1..hmm.m {
        fill_interior_transitions(cm, &map, &mut hmm, k, &psi);
    }
    fill_special_transitions(cm, &map, &mut hmm, &psi);
    hmm
}

/// The CP9's **own** expected occupancy `phi[k][match|insert|delete]` — a port of `fill_phi_cp9`
/// (start position `spos = 1`). A forward recursion through the built CP9: node 0's match (the
/// begin state B) is entered once; each later state accumulates occupancy from the previous
/// node's states via the CP9 transitions (plus `begin[k]`), and inserts fold their self-loop.
/// If the build is correct this reproduces the CM occupancy `psi` summed over the mapped states.
pub(crate) fn fill_phi(hmm: &Cp9) -> Vec<[f64; 3]> {
    let m = hmm.m;
    let mut phi = vec![[0.0f64; 3]; m + 1];
    let t = &hmm.t;

    phi[0][HMM_MATCH] = 1.0; // B (M_0) is entered in every parse
    phi[0][HMM_INSERT] = phi[0][HMM_MATCH] * t[0][CTMI];
    phi[0][HMM_INSERT] += phi[0][HMM_INSERT] * (t[0][CTII] / (1.0 - t[0][CTII]));
    phi[0][HMM_DELETE] = 0.0;

    for k in 1..=m {
        phi[k][HMM_MATCH] = phi[k - 1][HMM_MATCH] * t[k - 1][CTMM]
            + phi[k - 1][HMM_DELETE] * t[k - 1][CTDM]
            + phi[k - 1][HMM_INSERT] * t[k - 1][CTIM]
            + hmm.begin[k];
        // Deletes before inserts (an insert can follow this node's delete).
        phi[k][HMM_DELETE] = phi[k - 1][HMM_MATCH] * t[k - 1][CTMD]
            + phi[k - 1][HMM_INSERT] * t[k - 1][CTID]
            + phi[k - 1][HMM_DELETE] * t[k - 1][CTDD];
        phi[k][HMM_INSERT] = phi[k][HMM_MATCH] * t[k][CTMI] + phi[k][HMM_DELETE] * t[k][CTDI];
        phi[k][HMM_INSERT] += phi[k][HMM_INSERT] * (t[k][CTII] / (1.0 - t[k][CTII]));
    }
    phi
}

/// True if some base pair is modelled across two *adjacent* consensus columns (one MATP node
/// owning columns `k` and `k+1`) — a port of `check_cm_adj_bp`. The CP9 can't cleanly mirror that
/// node's single insert, so Infernal relaxes the psi-vs-phi check for such models.
pub(crate) fn has_adjacent_basepair(map: &Cp9Map) -> bool {
    (2..=map.hmm_m).any(|k| map.pos2nd[k] == map.pos2nd[k - 1])
}

/// The CP9 in **log-probability** form, ready for the Forward/Backward DP (`CP9Logoddsify`).
/// Emissions are log-probabilities rather than HMMER's log-*odds*; the null model cancels in the
/// state posteriors used for banding, so this changes the Forward score by a per-residue constant
/// but leaves the posteriors (hence the bands) identical, and keeps the DP self-contained.
pub(crate) struct Cp9Scores {
    pub m: usize,
    pub k: usize,
    pub tsc: Vec<[f64; 10]>, // [k][CT*]
    pub bsc: Vec<f64>,       // [k] begin (B→M_k)
    pub esc: Vec<f64>,       // [k] end (M_k→E)
    pub msc: Vec<Vec<f64>>,  // [k][residue] match emission
    pub isc: Vec<Vec<f64>>,  // [k][residue] insert emission
}

impl Cp9Scores {
    pub fn from_cp9(hmm: &Cp9) -> Cp9Scores {
        let m = hmm.m;
        let tsc = (0..=m)
            .map(|k| {
                let mut a = [f64::NEG_INFINITY; 10];
                for (c, slot) in a.iter_mut().enumerate() {
                    *slot = ln(hmm.t[k][c]);
                }
                a
            })
            .collect();
        let bsc = (0..=m).map(|k| ln(hmm.begin[k])).collect();
        let esc = (0..=m).map(|k| ln(hmm.end[k])).collect();
        let msc = (0..=m)
            .map(|k| (0..hmm.k).map(|x| ln(hmm.mat[k][x])).collect())
            .collect();
        let isc = (0..=m)
            .map(|k| (0..hmm.k).map(|x| ln(hmm.ins[k][x])).collect())
            .collect();
        Cp9Scores {
            m,
            k: hmm.k,
            tsc,
            bsc,
            esc,
            msc,
            isc,
        }
    }
}

/// A CP9 DP matrix: log-space scores for the match / insert / delete state of each node, at each
/// sequence row `0..=L`.
pub(crate) struct Cp9Mx {
    pub mmx: Vec<Vec<f64>>,
    pub imx: Vec<Vec<f64>>,
    pub dmx: Vec<Vec<f64>>,
}

/// CP9 **Forward** DP (glocal, full retained matrix) — a port of `cp9_Forward`. `dsq` is 1-based
/// with sentinels at `0` and `L+1`. Returns the matrix and the total log-probability of the
/// sequence (the score at the end state after the last residue).
pub(crate) fn cp9_forward(s: &Cp9Scores, dsq: &[Dsq]) -> (Cp9Mx, f64) {
    let m = s.m;
    let l = dsq.len() - 2;
    let ninf = f64::NEG_INFINITY;
    let mut mmx = vec![vec![ninf; m + 1]; l + 1];
    let mut imx = vec![vec![ninf; m + 1]; l + 1];
    let mut dmx = vec![vec![ninf; m + 1]; l + 1];

    // Row 0: B is entered (M_0 = 0); the all-delete path is reachable for every node.
    mmx[0][0] = 0.0;
    for k in 1..=m {
        dmx[0][k] = logsum3(
            mmx[0][k - 1] + s.tsc[k - 1][CTMD],
            imx[0][k - 1] + s.tsc[k - 1][CTID],
            dmx[0][k - 1] + s.tsc[k - 1][CTDD],
        );
    }

    for j in 1..=l {
        let res = dsq[j] as usize;
        imx[j][0] = s.isc[0][res]
            + logsum3(
                mmx[j - 1][0] + s.tsc[0][CTMI],
                imx[j - 1][0] + s.tsc[0][CTII],
                dmx[j - 1][0] + s.tsc[0][CTDI],
            );
        for k in 1..=m {
            mmx[j][k] = s.msc[k][res]
                + logsum4(
                    mmx[j - 1][k - 1] + s.tsc[k - 1][CTMM],
                    imx[j - 1][k - 1] + s.tsc[k - 1][CTIM],
                    dmx[j - 1][k - 1] + s.tsc[k - 1][CTDM],
                    mmx[j - 1][0] + s.bsc[k], // begin: B→M_k (only at the first emitting row)
                );
            imx[j][k] = s.isc[k][res]
                + logsum3(
                    mmx[j - 1][k] + s.tsc[k][CTMI],
                    imx[j - 1][k] + s.tsc[k][CTII],
                    dmx[j - 1][k] + s.tsc[k][CTDI],
                );
            dmx[j][k] = logsum3(
                mmx[j][k - 1] + s.tsc[k - 1][CTMD],
                imx[j][k - 1] + s.tsc[k - 1][CTID],
                dmx[j][k - 1] + s.tsc[k - 1][CTDD],
            );
        }
    }

    // End: any M_k→E (only M_M in glocal), plus D_M→E and I_M→E.
    let mut endsc = ninf;
    for k in 1..=m {
        endsc = logsum(endsc, mmx[l][k] + s.esc[k]);
    }
    endsc = logsum(endsc, dmx[l][m] + s.tsc[m][CTDM]);
    endsc = logsum(endsc, imx[l][m] + s.tsc[m][CTIM]);
    (Cp9Mx { mmx, imx, dmx }, endsc)
}

/// CP9 **Backward** DP (glocal, full retained matrix) — a port of `cp9_Backward`. Returns the
/// matrix (rows `1..=L` are the backward scores at emitting positions) and the total
/// log-probability (the score at the B state before the first residue), which must equal the
/// Forward score.
pub(crate) fn cp9_backward(s: &Cp9Scores, dsq: &[Dsq]) -> (Cp9Mx, f64) {
    let m = s.m;
    let l = dsq.len() - 2;
    let ninf = f64::NEG_INFINITY;
    let mut mmx = vec![vec![ninf; m + 1]; l + 1];
    let mut imx = vec![vec![ninf; m + 1]; l + 1];
    let mut dmx = vec![vec![ninf; m + 1]; l + 1];

    // Row L (last residue): everything ends in E.
    let rl = dsq[l] as usize;
    mmx[l][m] = s.esc[m] + s.msc[m][rl];
    imx[l][m] = s.tsc[m][CTIM] + s.isc[m][rl];
    dmx[l][m] = s.tsc[m][CTDM];
    for k in (1..m).rev() {
        mmx[l][k] = logsum(s.esc[k], dmx[l][k + 1] + s.tsc[k][CTMD]) + s.msc[k][rl];
        imx[l][k] = dmx[l][k + 1] + s.tsc[k][CTID] + s.isc[k][rl];
        dmx[l][k] = dmx[l][k + 1] + s.tsc[k][CTDD];
    }
    mmx[l][0] = dmx[l][1] + s.tsc[0][CTMD];
    imx[l][0] = dmx[l][1] + s.tsc[0][CTID] + s.isc[0][rl];

    // Rows L-1 .. 1.
    for i in (1..l).rev() {
        let res = dsq[i] as usize;
        mmx[i][m] = imx[i + 1][m] + s.tsc[m][CTMI] + s.msc[m][res];
        imx[i][m] = imx[i + 1][m] + s.tsc[m][CTII] + s.isc[m][res];
        dmx[i][m] = imx[i + 1][m] + s.tsc[m][CTDI];
        for k in (1..m).rev() {
            mmx[i][k] = logsum3(
                mmx[i + 1][k + 1] + s.tsc[k][CTMM],
                imx[i + 1][k] + s.tsc[k][CTMI],
                dmx[i][k + 1] + s.tsc[k][CTMD],
            ) + s.msc[k][res];
            imx[i][k] = logsum3(
                mmx[i + 1][k + 1] + s.tsc[k][CTIM],
                imx[i + 1][k] + s.tsc[k][CTII],
                dmx[i][k + 1] + s.tsc[k][CTID],
            ) + s.isc[k][res];
            dmx[i][k] = logsum3(
                mmx[i + 1][k + 1] + s.tsc[k][CTDM],
                imx[i + 1][k] + s.tsc[k][CTDI],
                dmx[i][k + 1] + s.tsc[k][CTDD],
            );
        }
        imx[i][0] = logsum3(
            mmx[i + 1][1] + s.tsc[0][CTIM],
            imx[i + 1][0] + s.tsc[0][CTII],
            dmx[i][1] + s.tsc[0][CTID],
        ) + s.isc[0][res];
        let mut b = ninf;
        for k in 1..=m {
            b = logsum(b, mmx[i + 1][k] + s.bsc[k]);
        }
        b = logsum(b, imx[i + 1][0] + s.tsc[0][CTMI]);
        b = logsum(b, dmx[i][1] + s.tsc[0][CTMD]);
        mmx[i][0] = b;
    }

    // Row 0 (B before the first residue): the total. Deletes here emit nothing; stored into the
    // matrix (dmx[0][k], mmx[0][0]) so the band computation can use the position-0 boundary.
    dmx[0][m] = imx[1][m] + s.tsc[m][CTDI];
    for k in (1..m).rev() {
        dmx[0][k] = logsum3(
            mmx[1][k + 1] + s.tsc[k][CTDM],
            imx[1][k] + s.tsc[k][CTDI],
            dmx[0][k + 1] + s.tsc[k][CTDD],
        );
    }
    let mut total = ninf;
    for k in 1..=m {
        total = logsum(total, mmx[1][k] + s.bsc[k]);
    }
    total = logsum(total, imx[1][0] + s.tsc[0][CTMI]);
    total = logsum(total, dmx[0][1] + s.tsc[0][CTMD]);
    mmx[0][0] = total; // B is occupied at position 0 in every parse
    (Cp9Mx { mmx, imx, dmx }, total)
}

/// CP9 state **posteriors** — `pmx[state][i][k]` = log P(node `k`'s state occupies position `i`
/// given the sequence), a port of `cp9_Posterior`. Combines Forward × Backward, dividing out the
/// emission (double-counted in both passes) and the total. Node 0's match is the non-emitting B
/// state, so its posterior is left at −∞. These feed the band computation (milestone 4).
pub(crate) fn cp9_posterior(
    s: &Cp9Scores,
    dsq: &[Dsq],
    fmx: &Cp9Mx,
    bmx: &Cp9Mx,
    total: f64,
) -> Cp9Mx {
    let m = s.m;
    let l = dsq.len() - 2;
    let ninf = f64::NEG_INFINITY;
    let mut mmx = vec![vec![ninf; m + 1]; l + 1];
    let mut imx = vec![vec![ninf; m + 1]; l + 1];
    let mut dmx = vec![vec![ninf; m + 1]; l + 1];
    for i in 1..=l {
        let res = dsq[i] as usize;
        for k in 1..=m {
            mmx[i][k] = fmx.mmx[i][k] + bmx.mmx[i][k] - s.msc[k][res] - total;
        }
        for k in 0..=m {
            imx[i][k] = fmx.imx[i][k] + bmx.imx[i][k] - s.isc[k][res] - total;
            dmx[i][k] = fmx.dmx[i][k] + bmx.dmx[i][k] - total;
        }
    }
    Cp9Mx { mmx, imx, dmx }
}

/// HMM **position bands** — for each CP9 node `k`, the band `[pn_min, pn_max]` of sequence
/// positions its match / insert / delete state can occupy (a port of `cp9_FB2HMMBands`). A state
/// with no posterior mass leaves both at `-1`.
pub(crate) struct HmmBands {
    pub m: usize,
    pub pn_min_m: Vec<i32>,
    pub pn_max_m: Vec<i32>,
    pub pn_min_i: Vec<i32>,
    pub pn_max_i: Vec<i32>,
    pub pn_min_d: Vec<i32>,
    pub pn_max_d: Vec<i32>,
}

/// Log posterior that node `k`'s match state occupies position `ip` (`0` = before the sequence,
/// where only the non-emitting B state lives).
fn post_m(
    s: &Cp9Scores,
    dsq: &[Dsq],
    f: &Cp9Mx,
    b: &Cp9Mx,
    total: f64,
    ip: usize,
    k: usize,
) -> f64 {
    if ip == 0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k == 0 {
        return f64::NEG_INFINITY;
    }
    f.mmx[ip][k] + b.mmx[ip][k] - s.msc[k][dsq[ip] as usize] - total
}

fn post_i(
    s: &Cp9Scores,
    dsq: &[Dsq],
    f: &Cp9Mx,
    b: &Cp9Mx,
    total: f64,
    ip: usize,
    k: usize,
) -> f64 {
    if ip == 0 {
        return f64::NEG_INFINITY; // inserts always emit
    }
    f.imx[ip][k] + b.imx[ip][k] - s.isc[k][dsq[ip] as usize] - total
}

fn post_d(
    _s: &Cp9Scores,
    _dsq: &[Dsq],
    f: &Cp9Mx,
    b: &Cp9Mx,
    total: f64,
    ip: usize,
    k: usize,
) -> f64 {
    f.dmx[ip][k] + b.dmx[ip][k] - total // deletes don't emit; row 0 is the all-delete prefix
}

/// Scan a state's per-position posterior for the band edges: `pn_min` = leftmost position whose
/// left-cumulative posterior mass exceeds the excluded-tail threshold, `pn_max` = rightmost from
/// the right. `-1` if the state never accumulates enough mass (effectively unused).
fn band_edges(l: usize, thresh: f64, post: impl Fn(usize) -> f64) -> (i32, i32) {
    let ninf = f64::NEG_INFINITY;
    let mut pn_min = -1;
    let mut mass = ninf;
    for ip in 0..=l {
        mass = logsum(mass, post(ip));
        if mass > thresh {
            pn_min = ip as i32;
            break;
        }
    }
    let mut pn_max = -1;
    mass = ninf;
    for ip in (0..=l).rev() {
        mass = logsum(mass, post(ip));
        if mass > thresh {
            pn_max = ip as i32;
            break;
        }
    }
    (pn_min, pn_max)
}

/// Posterior → HMM position bands (`cp9_FB2HMMBands`). `tau` is the per-side tail probability
/// excluded from each band (Infernal's `cm->tau`, default ~1e-7 → very wide, safe bands): the
/// band covers all but `tau` of each state's posterior mass.
pub(crate) fn cp9_fb2hmmbands(
    s: &Cp9Scores,
    dsq: &[Dsq],
    fmx: &Cp9Mx,
    bmx: &Cp9Mx,
    total: f64,
    tau: f64,
) -> HmmBands {
    let m = s.m;
    let l = dsq.len() - 2;
    let thresh = (tau / 2.0).ln(); // excluded mass on each side = (1 − (1−tau))/2 = tau/2
    let mut bands = HmmBands {
        m,
        pn_min_m: vec![-1; m + 1],
        pn_max_m: vec![-1; m + 1],
        pn_min_i: vec![-1; m + 1],
        pn_max_i: vec![-1; m + 1],
        pn_min_d: vec![-1; m + 1],
        pn_max_d: vec![-1; m + 1],
    };
    for k in 0..=m {
        let (mn, mx) = band_edges(l, thresh, |ip| post_m(s, dsq, fmx, bmx, total, ip, k));
        bands.pn_min_m[k] = mn;
        bands.pn_max_m[k] = mx;
        let (mn, mx) = band_edges(l, thresh, |ip| post_i(s, dsq, fmx, bmx, total, ip, k));
        bands.pn_min_i[k] = mn;
        bands.pn_max_i[k] = mx;
        let (mn, mx) = band_edges(l, thresh, |ip| post_d(s, dsq, fmx, bmx, total, ip, k));
        bands.pn_min_d[k] = mn;
        bands.pn_max_d[k] = mx;
    }
    bands
}

// ---------------------------------------------------------------------------
// Milestone 4b: HMM position bands -> CM per-state i,j bands -> per-j d-bands.
//
// Ports `HMMBandsEnforceValidParse` + `cp9_HMM2ijBands` + `ij2d_bands` from
// `hmmband.c`, restricted to the glocal case we run (no local begins/ends, no
// EL, `do_trunc = FALSE`): all the local/truncation/EL branches drop out.
// ---------------------------------------------------------------------------

/// Per-CM-state position bands derived from the CP9 HMM bands. `imin[v]..imax[v]` is the band on
/// the left border `i`, `jmin[v]..jmax[v]` the band on the right border `j` of the subsequence
/// emitted by the subtree rooted at state `v` (1-based positions). The sentinel for an unreachable
/// state is `imin=-1, imax=-2` (likewise for `j`). `hdmin[v][j-jmin[v]] .. hdmax[v][..]` is the
/// per-`(v,j)` band on the subsequence length `d` — the "doubly-ragged" structure the banded CYK
/// (milestone 5) consumes.
pub(crate) struct Cp9Bands {
    pub m: usize,
    pub imin: Vec<i32>,
    pub imax: Vec<i32>,
    pub jmin: Vec<i32>,
    pub jmax: Vec<i32>,
    pub hdmin: Vec<Vec<i32>>,
    pub hdmax: Vec<Vec<i32>>,
}

/// The "reachable residue" bands per HMM state, filled by [`hmm_bands_enforce_valid_parse`]. For
/// each HMM state, `r_*n[k]..r_*x[k]` is the set of residue positions for which that state is
/// reachable under *some* HMM parse that stays within the position bands. The `r_nn_i`/`r_nx_i`
/// (used to set CM `i` bands) and `r_nn_j`/`r_nx_j` (used for `j` bands) differ by an off-by-one on
/// the non-emitting delete/`M_0` states (see the long comment in the C source).
struct RBands {
    r_mn: Vec<i32>,
    r_mx: Vec<i32>,
    r_in: Vec<i32>,
    r_ix: Vec<i32>,
    r_dn: Vec<i32>,
    r_dx: Vec<i32>,
    r_nn_i: Vec<i32>,
    r_nx_i: Vec<i32>,
    r_nn_j: Vec<i32>,
    r_nx_j: Vec<i32>,
}

/// `HMMBandsFillGap`: two transitions into the same HMM state imply reachable ranges with a hole
/// between them; widen `I_k`'s position band so the hole can be bridged, letting us treat the
/// merged range as reachable. Mutates `pn_*_i[k]`. (`max2` is unused, mirroring the C signature.)
#[allow(clippy::too_many_arguments)]
fn hmm_bands_fill_gap(
    bands: &mut HmmBands,
    k: usize,
    min1: i32,
    min2: i32,
    prv_nd_r_mn: i32,
    prv_nd_r_dn: i32,
) {
    let right_min = if min1 <= min2 { min2 } else { min1 };
    let mut in_ = i32::MAX;
    if prv_nd_r_mn != i32::MAX {
        in_ = in_.min(prv_nd_r_mn + 1);
    }
    if prv_nd_r_dn != i32::MAX {
        in_ = in_.min(prv_nd_r_dn);
    }
    let ix = right_min - 1;
    bands.pn_min_i[k] = if bands.pn_min_i[k] != -1 {
        bands.pn_min_i[k].min(in_)
    } else {
        in_
    };
    bands.pn_max_i[k] = if bands.pn_max_i[k] != -1 {
        bands.pn_max_i[k].max(ix)
    } else {
        ix
    };
}

/// `HMMBandsFixUnreachable`: node `k` can't be reached within the current position bands. Greedily
/// widen the bands (scenario 1: let `I_{k-1}` emit the missing residues; scenario 2: force the
/// delete states of nodes `k..` to be used) so that at least one parse reaches `k`. Mutates the
/// `pn_*` bands. Only fires for impractically tight `tau`.
fn hmm_bands_fix_unreachable(bands: &mut HmmBands, k: usize, r_prv_min: i32) {
    let mut nxt_m = -1i32;
    let mut nxt_d = -1i32;
    if bands.pn_min_m[k] != -1 && bands.pn_max_m[k] != -1 && bands.pn_max_m[k] - 1 > r_prv_min {
        nxt_m = (bands.pn_min_m[k] - 1).max(r_prv_min);
    }
    if bands.pn_min_d[k] != -1 && bands.pn_max_d[k] != -1 && bands.pn_max_d[k] > r_prv_min {
        nxt_d = bands.pn_min_d[k].max(r_prv_min);
    }
    if nxt_m != -1 || nxt_d != -1 {
        // scenario 1: bridge the hole with I_{k-1}.
        let nxt_n = if nxt_m == -1 {
            nxt_d
        } else if nxt_d == -1 {
            nxt_m
        } else {
            nxt_m.min(nxt_d)
        };
        bands.pn_min_i[k - 1] = if bands.pn_min_i[k - 1] != -1 {
            bands.pn_min_i[k - 1].min(r_prv_min + 1)
        } else {
            r_prv_min + 1
        };
        bands.pn_max_i[k - 1] = if bands.pn_max_i[k - 1] != -1 {
            bands.pn_max_i[k - 1].max(nxt_n)
        } else {
            nxt_n
        };
    } else {
        // scenario 2: force delete states k..kp to absorb the over-emitted residues.
        let mut kp = k;
        while kp <= bands.m
            && bands.pn_max_m[kp] < (r_prv_min + 1)
            && bands.pn_max_d[kp] < r_prv_min
        {
            bands.pn_min_d[kp] = r_prv_min;
            bands.pn_max_d[kp] = r_prv_min;
            kp += 1;
        }
    }
}

/// Port of `HMMBandsEnforceValidParse` (glocal). Walks the HMM left-to-right propagating, for each
/// state, the residue positions reachable under the position bands; where the bands admit no parse
/// it widens them (via [`hmm_bands_fill_gap`] / [`hmm_bands_fix_unreachable`]) so a parse exists.
/// Returns the `r_*` reachable bands consumed by [`cp9_hmm2ij_bands`]. Mutates `bands` in place.
fn hmm_bands_enforce_valid_parse(
    s: &Cp9Scores,
    bands: &mut HmmBands,
    i0: i32,
    j0: i32,
    doing_search: bool,
) -> RBands {
    let hmm_m = bands.m;
    let big = i32::MAX;
    let small = i32::MIN;
    let n1 = hmm_m + 1;
    let mut r_mn = vec![big; n1];
    let mut r_mx = vec![small; n1];
    let mut r_in = vec![big; n1];
    let mut r_ix = vec![small; n1];
    let mut r_dn = vec![big; n1];
    let mut r_dx = vec![small; n1];
    let mut r_nn_i = vec![big; n1];
    let mut r_nx_i = vec![small; n1];
    let mut r_nn_j = vec![big; n1];
    let mut r_nx_j = vec![small; n1];
    let mut r_nn_hmm = vec![big; n1];
    let mut r_nx_hmm = vec![small; n1];
    let mut r_begn = big;
    let mut r_begx = small;
    let mut r_endn = big;
    let mut r_endx = small;
    let mut was_unr = vec![false; n1];
    let mut filled_gap = vec![false; n1];

    if bands.pn_min_m[0] != -1 {
        r_mn[0] = bands.pn_min_m[0];
        r_mx[0] = bands.pn_max_m[0];
    }

    // Transition relaxation into a target state's reachable band: project the source band forward by
    // `sd`, clip to the target's position band, and fold it in. Returns the (n, x) actually folded
    // (or None if there's no overlap). Inlined as a closure-free helper to keep the borrow simple.
    let mut k: i32 = 0;
    while k <= hmm_m as i32 {
        let ku = k as usize;
        let mut just_filled_gap = false;

        // --- transitions to the insert of node k (I_k) ---
        if r_mn[ku] <= r_mx[ku] && bands.pn_min_i[ku] != -1 {
            let n = r_mn[ku] + 1;
            let x = r_mx[ku] + 1;
            if x.min(bands.pn_max_i[ku]) - n.max(bands.pn_min_i[ku]) >= 0 {
                r_in[ku] = r_in[ku].min(n.max(bands.pn_min_i[ku]).min(bands.pn_max_i[ku]));
                r_ix[ku] = r_ix[ku].max(x.min(bands.pn_max_i[ku]));
            }
        }
        if r_dn[ku] <= r_dx[ku] && bands.pn_min_i[ku] != -1 {
            let n = r_dn[ku] + 1;
            let x = r_dx[ku] + 1;
            if x.min(bands.pn_max_i[ku]) - n.max(bands.pn_min_i[ku]) >= 0 {
                r_in[ku] = r_in[ku].min(n.max(bands.pn_min_i[ku]).min(bands.pn_max_i[ku]));
                r_ix[ku] = r_ix[ku].max(x.min(bands.pn_max_i[ku]));
            }
        }
        if r_in[ku] <= r_ix[ku] && bands.pn_min_i[ku] != -1 {
            // self-emitter: once entered we can emit up to pn_max_i[k].
            if r_in[ku] <= bands.pn_max_i[ku] {
                r_ix[ku] = bands.pn_max_i[ku];
            } else {
                r_in[ku] = big;
                r_ix[ku] = small;
            }
        }

        // --- transitions to the match of node k+1 (M_{k+1}) ---
        if k < hmm_m as i32 && bands.pn_min_m[ku + 1] != -1 {
            // sources: M_k, D_k, I_k (all project by +1 onto the match band).
            let srcs = [
                (r_mn[ku], r_mx[ku]),
                (r_dn[ku], r_dx[ku]),
                (r_in[ku], r_ix[ku]),
            ];
            for (sn, sx) in srcs {
                if sn > sx {
                    continue;
                }
                let n = sn + 1;
                let x = sx + 1;
                if x.min(bands.pn_max_m[ku + 1]) - n.max(bands.pn_min_m[ku + 1]) >= 0 {
                    let n = n.max(bands.pn_min_m[ku + 1]).min(bands.pn_max_m[ku + 1]);
                    let x = x.min(bands.pn_max_m[ku + 1]);
                    if r_mn[ku + 1] != big && x.min(r_mx[ku + 1]) - n.max(r_mn[ku + 1]) < -1 {
                        hmm_bands_fill_gap(bands, ku, n, r_mn[ku + 1], r_mn[ku - 1], r_dn[ku - 1]);
                        just_filled_gap = true;
                    }
                    r_mn[ku + 1] = r_mn[ku + 1].min(n);
                    r_mx[ku + 1] = r_mx[ku + 1].max(x);
                }
            }
        }
        // --- transitions to the END state (glocal: only from the last node) ---
        if k == hmm_m as i32 {
            if r_mn[ku] <= r_mx[ku] && s.esc[ku].is_finite() {
                r_endn = r_endn.min(r_mn[ku]);
                r_endx = r_endx.max(r_mx[ku]);
            }
            if r_dn[ku] <= r_dx[ku] && s.tsc[ku][CTDM].is_finite() {
                r_endn = r_endn.min(r_dn[ku]);
                r_endx = r_endx.max(r_dx[ku]);
            }
            if r_in[ku] <= r_ix[ku] && s.tsc[ku][CTIM].is_finite() {
                r_endn = r_endn.min(r_in[ku]);
                r_endx = r_endx.max(r_in[ku]); // C uses r_in[k] for both ends here
            }
        }

        // --- transitions to the delete of node k+1 (D_{k+1}) ---
        if k < hmm_m as i32 && bands.pn_min_d[ku + 1] != -1 {
            // sources: M_k, I_k, D_k (deletes don't emit, so no +1).
            let srcs = [
                (r_mn[ku], r_mx[ku]),
                (r_in[ku], r_ix[ku]),
                (r_dn[ku], r_dx[ku]),
            ];
            for (sn, sx) in srcs {
                if sn > sx {
                    continue;
                }
                let n = sn;
                let x = sx;
                if x.min(bands.pn_max_d[ku + 1]) - n.max(bands.pn_min_d[ku + 1]) >= 0 {
                    let n = n.max(bands.pn_min_d[ku + 1]).min(bands.pn_max_d[ku + 1]);
                    let x = x.min(bands.pn_max_d[ku + 1]);
                    if r_dn[ku + 1] != big && x.min(r_dx[ku + 1]) - n.max(r_dn[ku + 1]) < -1 {
                        hmm_bands_fill_gap(bands, ku, n, r_dn[ku + 1], r_mn[ku - 1], r_dn[ku - 1]);
                        just_filled_gap = true;
                    }
                    r_dn[ku + 1] = r_dn[ku + 1].min(n);
                    r_dx[ku + 1] = r_dx[ku + 1].max(x);
                }
            }
        }

        // --- update the reachable-by-node bands (r_nn_*, r_nx_*) ---
        if r_mn[ku] <= r_mx[ku] {
            r_nn_hmm[ku] = r_nn_hmm[ku].min(r_mn[ku]);
            r_nx_hmm[ku] = r_nx_hmm[ku].max(r_mx[ku]);
            if ku != hmm_m {
                r_nn_i[ku + 1] = r_nn_i[ku + 1].min(r_mn[ku] + 1);
                r_nx_i[ku + 1] = r_nx_i[ku + 1].max(r_mx[ku] + 1);
            }
            if ku != 0 {
                r_nn_j[ku - 1] = r_nn_j[ku - 1].min(r_mn[ku] - 1);
                r_nx_j[ku - 1] = r_nx_j[ku - 1].max(r_mx[ku] - 1);
            }
            if ku == hmm_m {
                if doing_search {
                    r_nn_j[ku] = r_nn_j[ku].min(r_mn[ku]);
                    r_nx_j[ku] = r_nx_j[ku].max(r_mx[ku]);
                } else if r_mx[ku] == j0 {
                    r_nn_j[ku] = r_nn_j[ku].min(j0);
                    r_nx_j[ku] = r_nx_j[ku].max(j0);
                }
            }
            if ku == 1 {
                if doing_search {
                    r_nn_i[ku] = r_nn_i[ku].min(r_mn[ku]);
                    r_nx_i[ku] = r_nx_i[ku].max(r_mx[ku]);
                    r_begn = r_begn.min(r_mn[ku]);
                    r_begx = r_begx.max(r_mx[ku]);
                } else if r_mn[ku] == i0 {
                    r_nn_i[ku] = r_nn_i[ku].min(i0);
                    r_nx_i[ku] = r_nx_i[ku].max(i0);
                }
            }
        }
        if r_in[ku] <= r_ix[ku] {
            r_nn_hmm[ku] = r_nn_hmm[ku].min(r_in[ku]);
            r_nx_hmm[ku] = r_nx_hmm[ku].max(r_ix[ku]);
            if ku != hmm_m {
                r_nn_i[ku + 1] = r_nn_i[ku + 1].min(r_in[ku] + 1);
                r_nx_i[ku + 1] = r_nx_i[ku + 1].max(r_ix[ku] + 1);
            }
            r_nn_j[ku] = r_nn_j[ku].min(r_in[ku] - 1);
            r_nx_j[ku] = r_nx_j[ku].max(r_ix[ku] - 1);
            if ku == 0 {
                r_begn = r_begn.min(r_in[ku]);
                r_begx = r_begx.max(r_ix[ku]);
            }
        }
        if r_dn[ku] <= r_dx[ku] {
            r_nn_hmm[ku] = r_nn_hmm[ku].min(r_dn[ku]);
            r_nx_hmm[ku] = r_nx_hmm[ku].max(r_dx[ku]);
            if ku != hmm_m {
                r_nn_i[ku + 1] = r_nn_i[ku + 1].min(r_dn[ku] + 1); // off-by-one
                r_nx_i[ku + 1] = r_nx_i[ku + 1].max(r_dx[ku] + 1);
            }
            if ku != 0 {
                r_nn_j[ku - 1] = r_nn_j[ku - 1].min(r_dn[ku]);
                r_nx_j[ku - 1] = r_nx_j[ku - 1].max(r_dx[ku]);
            }
            if ku == 1 {
                r_begn = r_begn.min(r_dn[ku] + 1);
                r_begx = r_begx.max(r_dx[ku] + 1);
            }
        }

        // --- is the node reachable? if not, widen bands and reprocess ---
        if r_mn[ku] > r_mx[ku] && r_dn[ku] > r_dx[ku] {
            assert!(ku != 0);
            assert!(!was_unr[ku], "cp9 bands: node {ku} unreachable on 2nd pass");
            was_unr[ku] = true;
            hmm_bands_fix_unreachable(bands, ku, r_nn_hmm[ku - 1]);
            k -= 2; // reprocess from k-1 after the loop's k += 1
        } else if just_filled_gap {
            assert!(
                !filled_gap[ku],
                "cp9 bands: node {ku} gap-filled on 2nd pass"
            );
            filled_gap[ku] = true;
            k -= 1; // reprocess k after the loop's k += 1
        }
        k += 1;
    }

    // Alignment pins the first/last emitted residue; search leaves the begin/end ranges as found.
    if !doing_search {
        r_begn = i0;
        r_begx = i0;
        r_endn = j0;
        r_endx = j0;
    }
    r_nn_i[1] = r_nn_i[1].min(r_begn);
    r_nx_i[1] = r_nx_i[1].max(r_begx);
    r_nn_j[hmm_m] = r_nn_j[hmm_m].min(r_endn);
    r_nx_j[hmm_m] = r_nx_j[hmm_m].max(r_endx);

    // Convert the INT_MAX/INT_MIN sentinels to the -1/-2 the band setter expects.
    let fix = |v: &mut [i32], lo: i32| {
        for x in v.iter_mut() {
            if *x == big {
                *x = -1;
            }
            if *x == small {
                *x = lo;
            }
        }
    };
    fix(&mut r_mn, -1);
    fix(&mut r_mx, -2);
    fix(&mut r_in, -1);
    fix(&mut r_ix, -2);
    fix(&mut r_dn, -1);
    fix(&mut r_dx, -2);
    fix(&mut r_nn_i, -1);
    fix(&mut r_nx_i, -2);
    fix(&mut r_nn_j, -1);
    fix(&mut r_nx_j, -2);

    RBands {
        r_mn,
        r_mx,
        r_in,
        r_ix,
        r_dn,
        r_dx,
        r_nn_i,
        r_nx_i,
        r_nn_j,
        r_nx_j,
    }
}

/// `ij2d_bands` (glocal, `do_trunc = FALSE`): from the per-state i,j bands, derive the per-`(v,j)`
/// band on the subsequence length `d`. `hdmin[v]`/`hdmax[v]` are indexed by `jp = j - jmin[v]`.
fn ij2d_bands(cm: &Cm, b: &mut Cp9Bands) {
    for v in 0..cm.m {
        let jspan = (b.jmax[v] - b.jmin[v]).max(-1); // -1 if unreachable -> empty
        let nj = (jspan + 1).max(0) as usize;
        let mut hdmin = vec![0i32; nj];
        let mut hdmax = vec![0i32; nj];
        if cm.sttype[v] == st::E {
            // E states always emit d=0.
            for jp in 0..nj {
                hdmin[jp] = 0;
                hdmax[jp] = 0;
            }
        } else {
            let sd = n_emit(cm.sttype[v]) as i32;
            for jp in 0..nj {
                let j = jp as i32 + b.jmin[v];
                let hdn = j - b.imax[v] + 1;
                let hdx = j - b.imin[v] + 1;
                if hdx < sd {
                    hdmin[jp] = -1;
                    hdmax[jp] = -2;
                } else {
                    hdmin[jp] = hdn.max(sd);
                    hdmax[jp] = hdx;
                }
            }
        }
        b.hdmin[v] = hdmin;
        b.hdmax[v] = hdmax;
    }
}

/// Port of `cp9_HMM2ijBands` (glocal, `do_trunc = FALSE`). Maps the HMM position bands onto CM
/// per-state i,j bands by walking the guide tree left-to-right (each node visited on the left then
/// the right, around the perimeter of the tree), then derives the per-j d-bands via [`ij2d_bands`].
/// `bands` is mutated by the enforce-valid-parse pass. `i0`/`j0` are the 1-based window endpoints.
pub(crate) fn cp9_hmm2ij_bands(
    cm: &Cm,
    map: &Cp9Map,
    s: &Cp9Scores,
    bands: &mut HmmBands,
    i0: i32,
    j0: i32,
    doing_search: bool,
) -> Cp9Bands {
    let hmm_m = bands.m;
    let r = hmm_bands_enforce_valid_parse(s, bands, i0, j0, doing_search);

    let m = cm.m;
    let mut imin = vec![-1i32; m];
    let mut imax = vec![-2i32; m];
    let mut jmin = vec![-1i32; m];
    let mut jmax = vec![-2i32; m];

    // Perimeter traversal: a stack of (on_right, node); a side-stack of remembered lpos for BIFs.
    let mut nd_pda: Vec<(bool, usize)> = vec![(false, 0)];
    let mut lpos_pda: Vec<i32> = Vec::new();
    let mut lpos: i32 = 0;
    let mut rpos: i32 = 0;

    while let Some((on_right, n)) = nd_pda.pop() {
        if on_right {
            match cm.ndtype[n] {
                nd::BIF => {
                    let v = cm.nodemap[n];
                    let w = cm.cfirst[v] as usize; // BEGL_S
                    let y = cm.cnum[v] as usize; // BEGR_S
                    imin[v] = if imin[w] != -1 { imin[w] } else { imin[y] };
                    imax[v] = if imax[w] != -2 { imax[w] } else { imax[y] };
                    jmin[v] = if jmin[y] != -1 { jmin[y] } else { jmin[w] };
                    jmax[v] = if jmax[y] != -2 { jmax[y] } else { jmax[w] };
                    // If either child is unreachable (only possible with local on, but kept), make
                    // the whole bifurcation unreachable.
                    if imin[w] == -1
                        || jmin[w] == -1
                        || imin[y] == -1
                        || jmin[y] == -1
                        || imax[w] == -2
                        || jmax[y] == -2
                    {
                        imin[v] = -1;
                        imin[w] = -1;
                        imin[y] = -1;
                        jmin[v] = -1;
                        jmin[w] = -1;
                        jmin[y] = -1;
                        imax[v] = -2;
                        imax[w] = -2;
                        imax[y] = -2;
                        jmax[v] = -2;
                        jmax[w] = -2;
                        jmax[y] = -2;
                        imin[y + 1] = -1;
                        jmin[y + 1] = -1;
                        imax[y + 1] = -2;
                        jmax[y + 1] = -2;
                    }
                }
                nd::MATP => {
                    lpos = map.nd2lpos[n];
                    rpos = map.nd2rpos[n];
                    let v = cm.nodemap[n];
                    let rp = rpos as usize;
                    jmin[v] = r.r_mn[rp]; // MATP_MP
                    jmax[v] = r.r_mx[rp];
                    jmin[v + 1] = r.r_dn[rp]; // MATP_ML
                    jmax[v + 1] = r.r_dx[rp];
                    jmin[v + 2] = r.r_mn[rp]; // MATP_MR
                    jmax[v + 2] = r.r_mx[rp];
                    jmin[v + 3] = r.r_dn[rp]; // MATP_D
                    jmax[v + 3] = r.r_dx[rp];
                    jmin[v + 4] = r.r_nn_j[rp - 1]; // MATP_IL
                    jmax[v + 4] = r.r_nx_j[rp - 1];
                    jmin[v + 5] = r.r_in[rp - 1]; // MATP_IR
                    jmax[v + 5] = r.r_ix[rp - 1];
                    imin[v + 5] = r.r_nn_i[(lpos + 1) as usize];
                    imax[v + 5] = r.r_nx_i[(lpos + 1) as usize];
                }
                nd::MATL => {
                    lpos = map.nd2lpos[n];
                    let v = cm.nodemap[n];
                    let rp = rpos as usize;
                    jmin[v] = r.r_nn_j[rp]; // MATL_ML
                    jmax[v] = r.r_nx_j[rp];
                    jmin[v + 1] = r.r_nn_j[rp]; // MATL_D
                    jmax[v + 1] = r.r_nx_j[rp];
                    jmin[v + 2] = r.r_nn_j[rp]; // MATL_IL
                    jmax[v + 2] = r.r_nx_j[rp];
                }
                nd::MATR => {
                    rpos = map.nd2rpos[n];
                    let v = cm.nodemap[n];
                    let rp = rpos as usize;
                    let lp = lpos as usize;
                    jmin[v] = r.r_mn[rp]; // MATR_MR
                    jmax[v] = r.r_mx[rp];
                    imin[v] = r.r_nn_i[lp];
                    imax[v] = r.r_nx_i[lp];
                    jmin[v + 1] = r.r_dn[rp]; // MATR_D
                    jmax[v + 1] = r.r_dx[rp];
                    imin[v + 1] = r.r_nn_i[lp];
                    imax[v + 1] = r.r_nx_i[lp];
                    jmin[v + 2] = r.r_in[rp - 1]; // MATR_IR
                    jmax[v + 2] = r.r_ix[rp - 1];
                    imin[v + 2] = r.r_nn_i[lp];
                    imax[v + 2] = r.r_nx_i[lp];
                }
                nd::BEGL | nd::BEGR => {
                    let v = cm.nodemap[n];
                    let mut imn = i32::MAX;
                    let mut imx = i32::MIN;
                    let mut jmn = i32::MAX;
                    let mut jmx = i32::MIN;
                    let cf = cm.cfirst[v] as usize;
                    let cn = cm.cnum[v] as usize;
                    for y in cf..cf + cn {
                        if imin[y] != -1 {
                            imn = imn.min(imin[y]);
                            imx = imx.max(imax[y]);
                        }
                        if jmin[y] != -1 {
                            jmn = jmn.min(jmin[y]);
                            jmx = jmx.max(jmax[y]);
                        }
                    }
                    if imn == i32::MAX {
                        imin[v] = -1;
                        imax[v] = -2;
                    } else {
                        imin[v] = imn;
                        imax[v] = imx;
                    }
                    if jmn == i32::MAX {
                        jmin[v] = -1;
                        jmax[v] = -2;
                    } else {
                        jmin[v] = jmn;
                        jmax[v] = jmx;
                    }
                    if cm.ndtype[n] == nd::BEGR {
                        // BEGR_IL: i band explicit from the HMM, emits before lpos.
                        let il = v + 1;
                        imin[il] = r.r_in[(lpos - 1) as usize];
                        imax[il] = r.r_ix[(lpos - 1) as usize];
                        if imin[il - 1] != -1 && imin[il] != -1 {
                            imin[il - 1] = imin[il - 1].min(imin[il]);
                            jmin[il] = if jmin[il - 1] == -1 {
                                -1
                            } else {
                                jmin[il - 1].max(i0)
                            };
                            jmax[il] = jmax[il - 1];
                        } else {
                            imin[il] = -1;
                            jmin[il] = -1;
                            imax[il] = -2;
                            jmax[il] = -2;
                        }
                        lpos = lpos_pda.pop().expect("lpos_pda underflow");
                    } else {
                        lpos_pda.push(lpos);
                        lpos = rpos + 1; // sister BEGR is on the right
                    }
                }
                nd::END => {
                    let v = cm.nodemap[n]; // END_E
                    let lp = lpos as usize;
                    imin[v] = r.r_nn_i[lp];
                    imax[v] = if r.r_nx_i[lp] == -2 {
                        -2
                    } else {
                        (r.r_nx_i[lp] + 1).min(j0 + 1) // +1 for StateDelta
                    };
                    if r.r_in[lp] != -1 {
                        imin[v] = imin[v].min((r.r_in[lp] - 1).max(i0));
                        imax[v] = imax[v].max((r.r_ix[lp] - 1).max(i0));
                    }
                    rpos = lpos;
                    if imin[v] != -1 {
                        jmin[v] = imin[v] - 1; // E emits d=0, so j = i-1
                        jmax[v] = imax[v] - 1;
                    }
                }
                nd::ROOT => {
                    debug_assert_eq!(lpos, 1);
                    debug_assert_eq!(rpos, hmm_m as i32);
                    let v = cm.nodemap[n]; // ROOT_S
                    imin[v] = r.r_nn_i[1];
                    imax[v] = r.r_nx_i[1];
                    jmin[v] = r.r_nn_j[hmm_m];
                    jmax[v] = r.r_nx_j[hmm_m];

                    let il = v + 1; // ROOT_IL
                    imin[il] = r.r_in[0];
                    imax[il] = r.r_ix[0];
                    jmin[il] = if r.r_nn_j[hmm_m] == -1 {
                        -1
                    } else {
                        r.r_nn_j[hmm_m].max(i0)
                    };
                    jmax[il] = r.r_nx_j[hmm_m];
                    if r.r_in[hmm_m] != -1 {
                        jmin[il] = jmin[il].min(r.r_in[hmm_m]);
                        jmax[il] = jmax[il].min(r.r_ix[hmm_m]);
                    }

                    let ir = v + 2; // ROOT_IR
                    if r.r_in[hmm_m] != -1 {
                        imin[ir] = r.r_nn_i[1];
                        imax[ir] = r.r_nx_i[1];
                        if imin[il] != -1 {
                            imin[ir] = imin[ir].min(imin[il] + 1);
                            imax[ir] = imax[ir].max(imax[il] + 1);
                        }
                        jmin[ir] = r.r_in[hmm_m];
                        jmax[ir] = r.r_ix[hmm_m];
                    }
                }
                _ => {}
            }
        } else {
            // on the left: set i bands for MATP / MATL only.
            match cm.ndtype[n] {
                nd::MATP => {
                    lpos = map.nd2lpos[n];
                    let v = cm.nodemap[n];
                    let lp = lpos as usize;
                    imin[v] = r.r_mn[lp]; // MATP_MP
                    imax[v] = r.r_mx[lp];
                    imin[v + 1] = r.r_mn[lp]; // MATP_ML
                    imax[v + 1] = r.r_mx[lp];
                    imin[v + 2] = if r.r_dn[lp] == -1 { -1 } else { r.r_dn[lp] + 1 }; // MATP_MR
                    imax[v + 2] = if r.r_dx[lp] == -2 { -2 } else { r.r_dx[lp] + 1 };
                    imin[v + 3] = if r.r_dn[lp] == -1 { -1 } else { r.r_dn[lp] + 1 }; // MATP_D
                    imax[v + 3] = if r.r_dx[lp] == -2 { -2 } else { r.r_dx[lp] + 1 };
                    imin[v + 4] = r.r_in[lp]; // MATP_IL
                    imax[v + 4] = r.r_ix[lp];
                    // MATP_IR i band is set on the right.
                }
                nd::MATL => {
                    lpos = map.nd2lpos[n];
                    let v = cm.nodemap[n];
                    let lp = lpos as usize;
                    imin[v] = r.r_mn[lp]; // MATL_ML
                    imax[v] = r.r_mx[lp];
                    imin[v + 1] = if r.r_dn[lp] == -1 { -1 } else { r.r_dn[lp] + 1 }; // MATL_D
                    imax[v + 1] = if r.r_dx[lp] == -2 { -2 } else { r.r_dx[lp] + 1 };
                    imin[v + 2] = r.r_in[lp]; // MATL_IL
                    imax[v + 2] = r.r_ix[lp];
                }
                _ => {}
            }

            // push children for the right pass (and recurse into subtree).
            if cm.ndtype[n] == nd::BIF {
                let v = cm.nodemap[n];
                let left_child = cm.ndidx[cm.cfirst[v] as usize];
                let right_child = cm.ndidx[cm.cnum[v] as usize];
                nd_pda.push((true, n)); // BIF right
                nd_pda.push((false, right_child)); // right child left
                nd_pda.push((false, left_child)); // left child left (popped first)
            } else {
                nd_pda.push((true, n)); // this node, right pass
                if cm.ndtype[n] != nd::END {
                    nd_pda.push((false, n + 1)); // child in tree order
                }
            }
        }
    }

    // Alignment: pin the root's i/j to span the whole window.
    if !doing_search {
        imin[0] = i0;
        if imin[1] != -1 {
            imin[1] = i0;
        }
        jmax[0] = j0;
        if jmin[1] != -1 {
            jmax[1] = j0;
        }
        if jmin[2] != -1 {
            jmax[2] = j0;
        }
    }

    // Final pass: collapse contradictory/detached states to unreachable; enforce minimum emission
    // for left emitters and MP states (glocal, do_trunc = FALSE).
    for v in 0..m {
        if imin[v] == -1 || jmin[v] == -1 {
            imin[v] = -1;
            jmin[v] = -1;
            imax[v] = -2;
            jmax[v] = -2;
        }
        if v + 1 < m && cm.stid[v + 1] == stid::END_E {
            // detached insert (an IL/MATP_IR immediately before an END)
            imin[v] = -1;
            jmin[v] = -1;
            imax[v] = -2;
            jmax[v] = -2;
        }
        if cm.sttype[v] == st::MP {
            if jmax[v] == i0 {
                imin[v] = -1;
                jmin[v] = -1;
                imax[v] = -2;
                jmax[v] = -2;
            } else if jmin[v] == i0 {
                jmin[v] += 1;
            }
        }
        if emits_left(cm.sttype[v]) && imin[v] != -1 {
            if jmax[v] == i0 - 1 {
                imin[v] = -1;
                jmin[v] = -1;
                imax[v] = -2;
                jmax[v] = -2;
            } else if jmin[v] == i0 - 1 {
                jmin[v] = i0;
            }
        }
    }

    let mut out = Cp9Bands {
        m,
        imin,
        imax,
        jmin,
        jmax,
        hdmin: vec![Vec::new(); m],
        hdmax: vec![Vec::new(); m],
    };
    ij2d_bands(cm, &mut out);
    out
}

/// Per-sequence driver: CP9 Forward/Backward → posteriors → HMM position bands → CM i,j,d bands.
/// `tau` is the per-side posterior tail excluded from each band; `doing_search` selects the search
/// vs. alignment band semantics. `dsq` is 1-based with sentinels at `0` and `L+1`.
pub(crate) fn cp9_seq2bands(
    cm: &Cm,
    map: &Cp9Map,
    s: &Cp9Scores,
    dsq: &[Dsq],
    tau: f64,
    doing_search: bool,
) -> Cp9Bands {
    let (fmx, ftot) = cp9_forward(s, dsq);
    let (bmx, _btot) = cp9_backward(s, dsq);
    let mut hmmbands = cp9_fb2hmmbands(s, dsq, &fmx, &bmx, ftot, tau);
    let l = (dsq.len() - 2) as i32;
    cp9_hmm2ij_bands(cm, map, s, &mut hmmbands, 1, l, doing_search)
}

// ---------------------------------------------------------------------------
// Milestone 5: doubly-ragged (j- and d-banded) CYK.
//
// A port of the HMM-banded CYK recursion from `cm_dpsearch.c::FastCYKScanHB`,
// in alignment mode (glocal, full-span ROOT->END over the window). The DP matrix
// `alpha[v][jp][dp]` allocates only in-band cells (`jp = j - jmin[v]`,
// `dp = d - hdmin[v][jp]`) — far smaller than the dense [v][j][d] grid. Given
// bands that contain the optimal parse (the milestone-4b no-drop guarantee), the
// root score is identical to the unbanded CYK.
// ---------------------------------------------------------------------------

#[inline]
fn cyk_add(a: f32, b: f32) -> f32 {
    if a <= IMPOSSIBLE || b <= IMPOSSIBLE {
        IMPOSSIBLE
    } else {
        a + b
    }
}

/// HMM-banded glocal CYK alignment score (bits) of the window `dsq` (1-based, sentinels at `0` and
/// `L+1`) under the i,j,d bands `b` (built in **alignment** mode for the same window). Returns the
/// root score `alpha[ROOT_S][j=L][d=L]`, byte-identical to the unbanded [`crate::align::align_glocal`]
/// when the bands contain the optimal parse. Returns `IMPOSSIBLE` if the full-span root cell is out
/// of band (no in-band parse spans the window).
pub(crate) fn cyk_banded_align(cm: &Cm, dsq: &[Dsq], b: &Cp9Bands) -> f32 {
    let l = (dsq.len() - 2) as i32;
    let m = cm.m;
    let oesc = Oesc::build(cm);
    let kp = oesc.kp;

    // Ragged DP matrix: alpha[v][jp][dp], only in-band cells. Unreachable states get no rows.
    let mut alpha: Vec<Vec<Vec<f32>>> = vec![Vec::new(); m];
    for v in 0..m {
        if b.jmin[v] < 0 {
            continue; // unreachable state
        }
        let njp = (b.jmax[v] - b.jmin[v] + 1) as usize;
        let mut rows = Vec::with_capacity(njp);
        for jp in 0..njp {
            let w = (b.hdmax[v][jp] - b.hdmin[v][jp] + 1).max(0) as usize;
            rows.push(vec![IMPOSSIBLE; w]);
        }
        alpha[v] = rows;
    }

    // In-band cell lookup: platonic alpha[y][j][d] -> ragged value, or None if out of band.
    let cell = |alpha: &[Vec<Vec<f32>>], y: usize, j: i32, d: i32| -> Option<f32> {
        if b.jmin[y] < 0 || j < b.jmin[y] || j > b.jmax[y] {
            return None;
        }
        let jp = (j - b.jmin[y]) as usize;
        if d < b.hdmin[y][jp] || d > b.hdmax[y][jp] {
            return None;
        }
        Some(alpha[y][jp][(d - b.hdmin[y][jp]) as usize])
    };

    for v in (0..m).rev() {
        if b.jmin[v] < 0 {
            continue;
        }
        let sttype = cm.sttype[v];
        let jminv = b.jmin[v];

        if sttype == st::E {
            for jp in 0..alpha[v].len() {
                // E: d == 0 only (hdmin == hdmax == 0).
                alpha[v][jp][0] = 0.0;
            }
            continue;
        }
        if sttype == st::B {
            let w = cm.cfirst[v] as usize; // left  (BEGL_S)
            let z = cm.cnum[v] as usize; // right (BEGR_S)
            if b.jmin[w] < 0 || b.jmin[z] < 0 {
                continue;
            }
            for jp_v in 0..alpha[v].len() {
                let j = jp_v as i32 + jminv;
                if j < b.jmin[z] || j > b.jmax[z] {
                    continue;
                }
                let hdminv = b.hdmin[v][jp_v];
                for dp_v in 0..alpha[v][jp_v].len() {
                    let d = dp_v as i32 + hdminv;
                    let mut best = IMPOSSIBLE;
                    // k = length of the right (z) fragment; left (w) spans j-k, d-k.
                    for k in 0..=d {
                        let (Some(left), Some(right)) =
                            (cell(&alpha, w, j - k, d - k), cell(&alpha, z, j, k))
                        else {
                            continue;
                        };
                        best = best.max(cyk_add(left, right));
                    }
                    alpha[v][jp_v][dp_v] = best;
                }
            }
            continue;
        }

        // ML, MR, MP, IL, IR, D, S: max over children (+ self for inserts), then emission.
        let sd = n_emit(sttype) as i32;
        let sdr = if emits_right(sttype) { 1 } else { 0 };
        let cfirst = cm.cfirst[v] as usize;
        let cnum = cm.cnum[v] as usize;
        let tsc = &cm.tsc[v];
        let single = oesc.single[v].as_deref();
        let pairtab = oesc.pair[v].as_deref();

        for jp_v in 0..alpha[v].len() {
            let j = jp_v as i32 + jminv;
            let hdminv = b.hdmin[v][jp_v];
            // d ascending so insert self-transitions read the finalized d-1 cell.
            for dp_v in 0..alpha[v][jp_v].len() {
                let d = dp_v as i32 + hdminv;
                let mut best = IMPOSSIBLE;
                for yoff in 0..cnum {
                    let y = cfirst + yoff;
                    if let Some(child) = cell(&alpha, y, j - sdr, d - sd) {
                        best = best.max(cyk_add(child, tsc[yoff]));
                    }
                }
                let i = j - d + 1; // window-local left position
                let emit = match sttype {
                    st::ML | st::IL => single.unwrap()[dsq[i as usize] as usize],
                    st::MR | st::IR => single.unwrap()[dsq[j as usize] as usize],
                    st::MP => {
                        pairtab.unwrap()[dsq[i as usize] as usize * kp + dsq[j as usize] as usize]
                    }
                    _ => 0.0, // D, S
                };
                alpha[v][jp_v][dp_v] = if sttype == st::D || sttype == st::S {
                    best
                } else {
                    cyk_add(best, emit)
                };
            }
        }
    }

    // Glocal alignment: root score is alpha[ROOT_S=0] at the full span (j = L, d = L).
    cell(&alpha, 0, l, l).unwrap_or(IMPOSSIBLE)
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
        let hmm = build_cp9(&cm);

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

    /// The CP9 build's correctness gate (`check_psi_vs_phi_cp9`): the CP9's own expected
    /// occupancy `phi[k][state]` must reproduce the CM occupancy `psi` summed over the CM states
    /// mapping to that CP9 node-state, within Infernal's 1e-4 threshold. This is the strong check
    /// that the subpath-summation transitions (and the emissions/map/begin/end) are all correct.
    #[test]
    fn cp9_psi_matches_phi() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);

        let map = build_cp9_map(&cm);
        let psi = expected_occupancy(&cm);
        let hmm = build_cp9(&cm);
        let phi = fill_phi(&hmm);

        // tRNA base pairs span loops, so no base pair is on adjacent columns → strict check.
        assert!(
            !has_adjacent_basepair(&map),
            "tRNA has no adjacent-column pairs"
        );

        let threshold = 1e-4;
        let mut violations = 0;
        for k in 0..=map.hmm_m {
            for (state, name) in [(HMM_MATCH, 'M'), (HMM_INSERT, 'I'), (HMM_DELETE, 'D')] {
                // Node 0 has no delete state in either model.
                let summed_psi: f64 = if k == 0 && state == HMM_DELETE {
                    0.0
                } else {
                    (0..2)
                        .filter_map(|s| {
                            let ap = map.hns2cs[k][state][s];
                            (ap >= 0).then(|| psi[ap as usize])
                        })
                        .sum()
                };
                let diff = (phi[k][state] - summed_psi).abs();
                if diff > threshold {
                    violations += 1;
                    eprintln!(
                        "{name} k={k}: phi={:.6} psi={:.6} diff={:.2e}",
                        phi[k][state], summed_psi, diff
                    );
                }
            }
        }
        assert_eq!(violations, 0, "{violations} psi-vs-phi violations (>1e-4)");
    }

    /// CP9 Forward/Backward consistency: the two DPs must agree on the total log-probability, and
    /// at every position the summed match+insert occupancy (`Σ_k F·B / emission`) must equal that
    /// total (a port of `cp9_CheckFB`). Also: every state posterior is a probability (≤ 1).
    #[test]
    fn cp9_forward_backward_consistent() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let scores = Cp9Scores::from_cp9(&build_cp9(&cm));

        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&cons).expect("digitize");
        let l = dsq.len() - 2;

        let (fmx, fsc) = cp9_forward(&scores, &dsq);
        let (bmx, bsc) = cp9_backward(&scores, &dsq);

        // 1. Forward and Backward agree on the total.
        assert!((fsc - bsc).abs() < 1e-6, "Forward {fsc} != Backward {bsc}");

        // 2. Per-position consistency + posteriors are probabilities.
        for i in 1..=l {
            let res = dsq[i] as usize;
            let mut fb = f64::NEG_INFINITY;
            // Match states emit at nodes 1..=M (node 0's match is the non-emitting B state).
            for k in 1..=scores.m {
                let pm = fmx.mmx[i][k] + bmx.mmx[i][k] - scores.msc[k][res];
                fb = logsum(fb, pm);
                assert!(
                    pm - fsc < 1e-6,
                    "match post i={i} k={k} = {}",
                    (pm - fsc).exp()
                );
            }
            // Insert states emit at nodes 0..=M (node 0's insert is ROOT_IL).
            for k in 0..=scores.m {
                let pi = fmx.imx[i][k] + bmx.imx[i][k] - scores.isc[k][res];
                fb = logsum(fb, pi);
                assert!(
                    pi - fsc < 1e-6,
                    "insert post i={i} k={k} = {}",
                    (pi - fsc).exp()
                );
            }
            assert!(
                (fb - fsc).abs() < 1e-4,
                "position {i}: Σ F·B = {fb} != total {fsc}"
            );
        }

        // 3. cp9_posterior: at each position the emitting (match k≥1 + insert) posteriors sum to 1.
        let pmx = cp9_posterior(&scores, &dsq, &fmx, &bmx, fsc);
        for i in 1..=l {
            let mut p = 0.0;
            for k in 1..=scores.m {
                p += pmx.mmx[i][k].exp();
            }
            for k in 0..=scores.m {
                p += pmx.imx[i][k].exp();
            }
            assert!(
                (p - 1.0).abs() < 1e-4,
                "position {i}: posterior mass {p} != 1"
            );
        }
    }

    /// HMM position bands on the model consensus: since the consensus aligns 1:1 to the model
    /// (length == M), node `k`'s match state emits position `k`, so the band `[pn_min_m, pn_max_m]`
    /// must contain `k`. Bands are also well-formed (`min ≤ max`).
    #[test]
    fn cp9_hmm_bands_contain_consensus_diagonal() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let scores = Cp9Scores::from_cp9(&build_cp9(&cm));

        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&cons).expect("digitize");
        let l = dsq.len() - 2;
        assert_eq!(
            l, scores.m,
            "consensus length should equal the model length M"
        );

        let (fmx, fsc) = cp9_forward(&scores, &dsq);
        let (bmx, _) = cp9_backward(&scores, &dsq);
        let bands = cp9_fb2hmmbands(&scores, &dsq, &fmx, &bmx, fsc, 1e-7);

        for k in 1..=scores.m {
            let (mn, mx) = (bands.pn_min_m[k], bands.pn_max_m[k]);
            assert!(mn >= 0 && mx >= 0, "node {k} match band unset");
            assert!(mn <= mx, "node {k} match band inverted {mn}..{mx}");
            assert!(
                mn <= k as i32 && k as i32 <= mx,
                "node {k} match band {mn}..{mx} excludes diagonal position {k}"
            );
        }
    }

    /// The no-drop guarantee (milestone 4b): the CM i,j,d bands derived from the HMM bands must
    /// **contain the true CYK parse**. We align the model consensus, then assert every parsed state
    /// `(v, i, j, d)` falls inside its band — if any didn't, HMM-banded search could silently lose
    /// the optimal hit.
    #[test]
    fn cp9_ij_bands_contain_consensus_parse() {
        use crate::align::align_glocal_parse;

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let map = build_cp9_map(&cm);
        let scores = Cp9Scores::from_cp9(&build_cp9(&cm));

        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&cons).expect("digitize");
        let l = dsq.len() - 2;

        let emap = EmitMap::build(&cm);
        let parse = align_glocal_parse(&cm, &dsq, 1, l, &emap);

        // Both band modes must contain the true full-span parse: alignment mode (the parse must
        // span the whole window) and search mode (the semantics milestone 5's scan will use).
        for doing_search in [false, true] {
            let bands = cp9_seq2bands(&cm, &map, &scores, &dsq, 1e-7, doing_search);
            for &(v, i, j, d) in &parse {
                let (i, j, d) = (i as i32, j as i32, d as i32);
                assert!(
                    bands.imin[v] != -1,
                    "[search={doing_search}] state {v} ({}) used by parse but marked unreachable",
                    cm.sttype[v]
                );
                assert!(
                    bands.imin[v] <= i && i <= bands.imax[v],
                    "[search={doing_search}] state {v} ({}): parse i={i} outside i-band {}..{}",
                    cm.sttype[v],
                    bands.imin[v],
                    bands.imax[v]
                );
                assert!(
                    bands.jmin[v] <= j && j <= bands.jmax[v],
                    "[search={doing_search}] state {v} ({}): parse j={j} outside j-band {}..{}",
                    cm.sttype[v],
                    bands.jmin[v],
                    bands.jmax[v]
                );
                let jp = (j - bands.jmin[v]) as usize;
                assert!(
                    bands.hdmin[v][jp] <= d && d <= bands.hdmax[v][jp],
                    "[search={doing_search}] state {v} ({}): parse d={d} outside d-band {}..{} at j={j}",
                    cm.sttype[v],
                    bands.hdmin[v][jp],
                    bands.hdmax[v][jp]
                );
            }
        }
    }

    /// No-drop on a *real* tRNA (not the consensus): its parse has indels relative to the model, so
    /// it pushes off the band diagonal and exercises the off-diagonal / reachability logic. We find
    /// the hit window with the unbanded CYK scanner, then check its parse against the HMM bands.
    #[test]
    fn cp9_ij_bands_contain_real_trna_parse() {
        use crate::align::align_glocal_parse;
        use crate::search::cyk_search;

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let map = build_cp9_map(&cm);
        let scores = Cp9Scores::from_cp9(&build_cp9(&cm));
        let emap = EmitMap::build(&cm);

        let seq = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_embedded.fa"
        ))
        .expect("read embedded")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&seq).expect("digitize");

        // Locate the tRNA with the unbanded glocal CYK scanner; take the best forward-strand hit.
        let hits = cyk_search(&cm, &dsq, cm.w as usize, 10.0);
        let hit = hits
            .iter()
            .filter(|h| matches!(h.strand, crate::search::Strand::Plus))
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
            .expect("a forward-strand tRNA hit");
        let (i0, j0) = (hit.i, hit.j);

        // Window sub-dsq (1-based, own sentinels) for the CP9 Forward/Backward over just the hit.
        let mut wdsq = vec![dsq[0]];
        wdsq.extend_from_slice(&dsq[i0..=j0]);
        wdsq.push(dsq[dsq.len() - 1]);

        // Parse is in window-local coordinates (1..=wlen); bands likewise, so they line up.
        let parse = align_glocal_parse(&cm, &dsq, i0, j0, &emap);
        let bands = cp9_seq2bands(&cm, &map, &scores, &wdsq, 1e-7, false);

        for &(v, i, j, d) in &parse {
            let (i, j, d) = (i as i32, j as i32, d as i32);
            assert!(
                bands.imin[v] != -1,
                "state {v} ({}) used by parse but marked unreachable",
                cm.sttype[v]
            );
            assert!(
                bands.imin[v] <= i && i <= bands.imax[v],
                "state {v} ({}): parse i={i} outside i-band {}..{}",
                cm.sttype[v],
                bands.imin[v],
                bands.imax[v]
            );
            assert!(
                bands.jmin[v] <= j && j <= bands.jmax[v],
                "state {v} ({}): parse j={j} outside j-band {}..{}",
                cm.sttype[v],
                bands.jmin[v],
                bands.jmax[v]
            );
            let jp = (j - bands.jmin[v]) as usize;
            assert!(
                bands.hdmin[v][jp] <= d && d <= bands.hdmax[v][jp],
                "state {v} ({}): parse d={d} outside d-band {}..{} at j={j}",
                cm.sttype[v],
                bands.hdmin[v][jp],
                bands.hdmax[v][jp]
            );
        }
    }

    /// Milestone 5: the doubly-ragged HMM-banded CYK must produce a score **byte-identical** to the
    /// unbanded CYK on both the consensus and a real tRNA — the bands (which contain the optimal
    /// parse, per 4b) prune the matrix without changing the optimum.
    #[test]
    fn cyk_banded_matches_unbanded() {
        use crate::align::align_glocal;
        use crate::search::cyk_search;

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let map = build_cp9_map(&cm);
        let scores = Cp9Scores::from_cp9(&build_cp9(&cm));
        let emap = EmitMap::build(&cm);

        // Build a window sub-dsq (1-based, own sentinels) and compare banded vs unbanded CYK on it.
        let check = |cm: &Cm, wdsq: &[Dsq]| {
            let l = wdsq.len() - 2;
            let unbanded = align_glocal(cm, wdsq, 1, l, &emap).score;
            let bands = cp9_seq2bands(cm, &map, &scores, wdsq, 1e-7, false);
            let banded = cyk_banded_align(cm, wdsq, &bands);
            assert_eq!(
                banded.to_bits(),
                unbanded.to_bits(),
                "banded CYK {banded} != unbanded {unbanded}"
            );
        };

        // 1. The model consensus (the whole sequence is the window).
        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        check(&cm, &cm.abc.digitize(&cons).expect("digitize"));

        // 2. A real embedded tRNA: locate it with the unbanded scanner, then window it out.
        let seq = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_embedded.fa"
        ))
        .expect("read embedded")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&seq).expect("digitize");
        let hit = cyk_search(&cm, &dsq, cm.w as usize, 10.0)
            .into_iter()
            .filter(|h| matches!(h.strand, crate::search::Strand::Plus))
            .max_by(|a, c| a.score.partial_cmp(&c.score).unwrap())
            .expect("a forward-strand tRNA hit");
        let mut wdsq = vec![dsq[0]];
        wdsq.extend_from_slice(&dsq[hit.i..=hit.j]);
        wdsq.push(dsq[dsq.len() - 1]);
        check(&cm, &wdsq);
    }
}
