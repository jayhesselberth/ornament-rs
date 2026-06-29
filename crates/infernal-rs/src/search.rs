//! Non-banded scanning DP (glocal), ported from `RefCYKScan` / `RefIInsideScan` in
//! `src/cm_dpsearch.c`.
//!
//! This is the correctness core: for every end position `j` it computes the best parse
//! score of every subsequence `i..j` rooted at ROOT_S, and reports the best-scoring hit.
//! The recursion is shared between **CYK** (Viterbi-SCFG, `⊕ = max`) and **Inside**
//! (`⊕ = log-sum-exp`) via the [`Semiring`] parameter — the only difference is how
//! alternative sub-parses are combined within a cell. Glocal mode (no local begins/ends)
//! is the first parity target, matching `cmsearch -g --max [--cyk]`.
//!
//! Memory is bounded as Infernal bounds it: two rolling `j` rows for ordinary states and
//! `W+1` rolling rows for `BEGL_S` states.

use easel_rs::Dsq;

use crate::config::IMPOSSIBLE;
use crate::emit::Oesc;
use crate::model::{emits_right, n_emit, st, stid, Cm};

/// A single hit (1-based, inclusive coordinates on the forward strand).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CykHit {
    pub score: f32,
    /// Start position (1-based).
    pub i: usize,
    /// End position (1-based).
    pub j: usize,
}

/// Result of a glocal scan.
#[derive(Debug, Clone)]
pub struct CykScan {
    /// Best hit found (None if the sequence is empty or nothing scored above IMPOSSIBLE).
    pub best: Option<CykHit>,
    /// Best ROOT_S score at each end position `j` (1-based; index 0 unused).
    pub bestsc_per_j: Vec<f32>,
}

/// The "⊕" operation that combines alternative sub-parses within a DP cell.
pub trait Semiring {
    /// Combine two scores (both in bits / log2 space).
    fn or(a: f32, b: f32) -> f32;
}

/// Max-plus semiring → CYK (best single parse).
pub struct MaxPlus;
impl Semiring for MaxPlus {
    #[inline]
    fn or(a: f32, b: f32) -> f32 {
        if a >= b {
            a
        } else {
            b
        }
    }
}

/// Log-sum-exp semiring (base 2) → Inside (sum over all parses).
pub struct LogSum;
impl Semiring for LogSum {
    #[inline]
    fn or(a: f32, b: f32) -> f32 {
        flogsum(a, b)
    }
}

/// `log2(2^a + 2^b)`, with the IMPOSSIBLE sentinel acting as `-inf`.
#[inline]
fn flogsum(a: f32, b: f32) -> f32 {
    if a <= IMPOSSIBLE {
        return b;
    }
    if b <= IMPOSSIBLE {
        return a;
    }
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    // hi + log2(1 + 2^(lo-hi)); computed in f64 for accuracy.
    hi + ((1.0f64 + ((lo - hi) as f64 * std::f64::consts::LN_2).exp()).ln()
        * std::f64::consts::LOG2_E) as f32
}

/// Saturating add (⊗): keeps the IMPOSSIBLE sentinel absorbing.
#[inline]
fn add(a: f32, b: f32) -> f32 {
    if a <= IMPOSSIBLE || b <= IMPOSSIBLE {
        IMPOSSIBLE
    } else {
        a + b
    }
}

/// Glocal CYK scan (best single parse). Coordinates of the best hit are reported.
pub fn cyk_scan_glocal(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan_glocal::<MaxPlus>(cm, dsq, w_max)
}

/// Glocal Inside scan (sum over all parses). The reported hit score is the Inside score
/// at the best end position; coordinates mark that endpoint.
pub fn inside_scan_glocal(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan_glocal::<LogSum>(cm, dsq, w_max)
}

/// Shared glocal scanning recursion, parameterized by the within-cell ⊕ ([`Semiring`]).
fn scan_glocal<S: Semiring>(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    let l = dsq.len().saturating_sub(2);
    let w = w_max.min(l);
    let mut bestsc_per_j = vec![IMPOSSIBLE; l + 1];
    if l == 0 {
        return CykScan { best: None, bestsc_per_j };
    }

    let oesc = Oesc::build(cm);
    let kp = oesc.kp;
    let m = cm.m;

    // Compact indexing for BEGL_S states (the only states needing the W+1 rolling deck).
    let mut begl_idx = vec![usize::MAX; m];
    let mut n_begl = 0;
    for v in 0..m {
        if cm.stid[v] == stid::BEGL_S {
            begl_idx[v] = n_begl;
            n_begl += 1;
        }
    }

    let width = w + 1;
    let mut alpha = vec![vec![IMPOSSIBLE; m * width]; 2];
    let mut alpha_begl = vec![vec![IMPOSSIBLE; n_begl.max(1) * width]; width];

    // --- d=0 base case (empty-fragment scores), per cm_scan_mx_InitializeFloats. ---
    // Only E/S/D/B states get a finite d=0; emitters stay IMPOSSIBLE. Glocal => endsc is
    // IMPOSSIBLE everywhere (no local ends).
    for v in (0..m).rev() {
        if cm.stid[v] == stid::BEGL_S {
            let ych = cm.cfirst[v] as usize;
            let mut s = IMPOSSIBLE; // endsc[v]
            for yoff in 0..cm.cnum[v] as usize {
                s = S::or(s, add(alpha[0][(ych + yoff) * width], cm.tsc[v][yoff]));
            }
            let bi = begl_idx[v];
            for jp in 0..width {
                alpha_begl[jp][bi * width] = s;
            }
        } else if cm.sttype[v] == st::E {
            alpha[0][v * width] = 0.0;
            alpha[1][v * width] = 0.0;
        } else if cm.sttype[v] == st::S || cm.sttype[v] == st::D {
            let ych = cm.cfirst[v] as usize;
            let mut s = IMPOSSIBLE;
            for yoff in 0..cm.cnum[v] as usize {
                s = S::or(s, add(alpha[0][(ych + yoff) * width], cm.tsc[v][yoff]));
            }
            alpha[0][v * width] = s;
            alpha[1][v * width] = s;
        } else if cm.sttype[v] == st::B {
            let wbi = begl_idx[cm.cfirst[v] as usize];
            let ych = cm.cnum[v] as usize;
            let s = add(alpha_begl[0][wbi * width], alpha[0][ych * width]);
            alpha[0][v * width] = s;
            alpha[1][v * width] = s;
        }
        // emitters: d=0 stays IMPOSSIBLE.
    }

    let mut best: Option<CykHit> = None;

    for j in 1..=l {
        let cur = j % 2;
        let prv = (j - 1) % 2;
        let dx = j.min(w);

        for v in (1..m).rev() {
            let sttype = cm.sttype[v];
            if sttype == st::E {
                continue; // constant, set in base case
            }
            let sd = n_emit(sttype);
            let right = emits_right(sttype);
            let is_begl = cm.stid[v] == stid::BEGL_S;

            if sttype == st::B {
                let wbi = begl_idx[cm.cfirst[v] as usize];
                let ych = cm.cnum[v] as usize;
                for d in 1..=dx {
                    let mut sc = IMPOSSIBLE;
                    for k in 0..=d {
                        let left = alpha_begl[(j - k) % width][wbi * width + (d - k)];
                        let rightsc = alpha[cur][ych * width + k];
                        sc = S::or(sc, add(left, rightsc));
                    }
                    alpha[cur][v * width + d] = sc;
                }
                continue;
            }

            let ych = cm.cfirst[v] as usize;
            let cnum = cm.cnum[v] as usize;
            let tsc = &cm.tsc[v];
            let child_row = if right { prv } else { cur };
            let single = oesc.single[v].as_deref();
            let pairtab = oesc.pair[v].as_deref();

            for d in 1..=dx {
                if d < sd {
                    continue; // emitters with too few residues stay IMPOSSIBLE
                }
                let mut sc = IMPOSSIBLE;
                let cd = d - sd;
                for yoff in 0..cnum {
                    sc = S::or(sc, add(alpha[child_row][(ych + yoff) * width + cd], tsc[yoff]));
                }
                let val = match sttype {
                    st::ML | st::IL => {
                        let i = j - d + 1;
                        add(sc, single.unwrap()[dsq[i] as usize])
                    }
                    st::MR | st::IR => add(sc, single.unwrap()[dsq[j] as usize]),
                    st::MP => {
                        let i = j - d + 1;
                        add(sc, pairtab.unwrap()[dsq[i] as usize * kp + dsq[j] as usize])
                    }
                    _ => sc, // D, S (EMITNONE)
                };
                if is_begl {
                    alpha_begl[j % width][begl_idx[v] * width + d] = val;
                } else {
                    alpha[cur][v * width + d] = val;
                }
            }
        }

        // ROOT_S (v = 0): glocal, transitions only. Report hits d >= 1.
        let tsc0 = &cm.tsc[0];
        let ych = cm.cfirst[0] as usize;
        let cnum0 = cm.cnum[0] as usize;
        let mut bestsc_j = IMPOSSIBLE;
        for d in 1..=dx {
            let mut sc = IMPOSSIBLE;
            for yoff in 0..cnum0 {
                sc = S::or(sc, add(alpha[cur][(ych + yoff) * width + d], tsc0[yoff]));
            }
            alpha[cur][d] = sc; // ROOT stored at v=0 offset (read by nobody)
            if sc > bestsc_j {
                bestsc_j = sc;
            }
            if best.map_or(true, |b| sc > b.score) {
                best = Some(CykHit {
                    score: sc,
                    i: j - d + 1,
                    j,
                });
            }
        }
        bestsc_per_j[j] = bestsc_j;
    }

    CykScan { best, bestsc_per_j }
}
