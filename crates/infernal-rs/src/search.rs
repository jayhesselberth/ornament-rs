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
use crate::evalues::{evalue, SearchMode};
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

/// Which strand a hit was found on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Forward (top) strand.
    Plus,
    /// Reverse-complement (bottom) strand.
    Minus,
}

/// A reported hit in ORIGINAL sequence coordinates, mirroring `CM_HIT` / `cmsearch --tblout`.
///
/// `i`/`j` follow the tblout `seq from`/`seq to` convention: on the [`Strand::Plus`] strand
/// `i < j`; on [`Strand::Minus`] the alignment runs 3'→5' along the top strand so `i > j`
/// (`i` is the high coordinate, `j` the low). The matched span is always `i.min(j)..=i.max(j)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Bit score.
    pub score: f32,
    /// E-value for the search mode, if the model is calibrated for it.
    pub evalue: Option<f64>,
    /// Start coordinate (1-based) on the original (forward) sequence.
    pub i: usize,
    /// End coordinate (1-based) on the original (forward) sequence.
    pub j: usize,
    /// Strand the hit was found on.
    pub strand: Strand,
}

/// A raw hit on the strand that was scanned, in that strand's local 1-based coordinates.
#[derive(Debug, Clone, Copy)]
struct RawHit {
    score: f32,
    i: usize,
    j: usize,
}

/// Gamma semi-HMM for optimal resolution of overlapping hits, ported from
/// `CreateGammaHitMx` / `UpdateGammaHitMx` / `TBackGammaHitMx` (`src/cm_mx.c`). This is the
/// non-greedy (default) Infernal path: it picks the maximum-total-score set of
/// non-overlapping hits, then a traceback recovers each above-cutoff hit.
struct Gamma {
    /// `mx[j]` = best cumulative score of all chosen hits ending at or before `j`.
    mx: Vec<f32>,
    /// `gback[j]` = start `i` of the hit ending at `j` (0 = no hit ends here).
    gback: Vec<usize>,
    /// `savesc[j]` = bit score of the hit ending at `j`.
    savesc: Vec<f32>,
    /// Reporting cutoff (bits); a recovered hit is kept iff `savesc >= cutoff`.
    cutoff: f32,
}

impl Gamma {
    fn new(l: usize, cutoff: f32) -> Self {
        Gamma {
            mx: vec![0.0; l + 1],
            gback: vec![0; l + 1],
            savesc: vec![IMPOSSIBLE; l + 1],
            cutoff,
        }
    }

    /// Traceback recovering all above-cutoff non-overlapping hits (`TBackGammaHitMx`),
    /// returned in ascending coordinate order.
    fn traceback(&self, l: usize) -> Vec<RawHit> {
        let mut hits = Vec::new();
        let mut j = l;
        while j >= 1 {
            let i = self.gback[j];
            if i == 0 {
                j -= 1; // no hit ends at j
            } else {
                if self.savesc[j] >= self.cutoff {
                    hits.push(RawHit {
                        score: self.savesc[j],
                        i,
                        j,
                    });
                }
                if i <= 1 {
                    break;
                }
                j = i - 1;
            }
        }
        hits.reverse();
        hits
    }
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

/// CYK scan (best single parse). Honours the model's configured mode: glocal if the CM was
/// configured with `configure_scores`, local if with `configure_local`.
pub fn cyk_scan(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan::<MaxPlus>(cm, dsq, w_max)
}

/// Inside scan (sum over all parses), mode per the CM's configuration.
pub fn inside_scan(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan::<LogSum>(cm, dsq, w_max)
}

/// Glocal CYK scan. Alias for [`cyk_scan`]; the CM must be configured glocal.
pub fn cyk_scan_glocal(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan::<MaxPlus>(cm, dsq, w_max)
}

/// Glocal Inside scan. Alias for [`inside_scan`]; the CM must be configured glocal.
pub fn inside_scan_glocal(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan::<LogSum>(cm, dsq, w_max)
}

/// CYK search of BOTH strands, returning all hits scoring `>= cutoff_bits`, with overlapping
/// windows resolved by the gamma semi-HMM (Infernal's default), sorted best-first (by
/// E-value when calibrated, else by score). Honours the CM's configured glocal/local mode.
pub fn cyk_search(cm: &Cm, dsq: &[Dsq], w_max: usize, cutoff_bits: f32) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalCyk
    } else {
        SearchMode::GlocalCyk
    };
    search::<MaxPlus>(cm, dsq, w_max, cutoff_bits, mode)
}

/// Inside counterpart of [`cyk_search`] (sum-over-parses scoring).
pub fn inside_search(cm: &Cm, dsq: &[Dsq], w_max: usize, cutoff_bits: f32) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalInside
    } else {
        SearchMode::GlocalInside
    };
    search::<LogSum>(cm, dsq, w_max, cutoff_bits, mode)
}

/// Both-strand multi-hit driver shared by [`cyk_search`]/[`inside_search`].
///
/// Scans the forward strand, then the reverse complement, mapping minus-strand windows back
/// to original coordinates. E-values use the total residues searched on both strands
/// (`2*L`), matching `cmsearch`'s default reporting.
fn search<S: Semiring>(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    mode: SearchMode,
) -> Vec<Hit> {
    let l = dsq.len().saturating_sub(2);
    let searched = 2.0 * l as f64;
    let mut hits: Vec<Hit> = Vec::new();

    // Forward strand.
    for h in scan_core::<S>(cm, dsq, w_max, Some(cutoff_bits)).1 {
        hits.push(Hit {
            score: h.score,
            evalue: evalue(cm, mode, h.score, searched),
            i: h.i,
            j: h.j,
            strand: Strand::Plus,
        });
    }

    // Reverse-complement strand. A window (a..b) in the revcomp (a<=b) maps to original
    // coordinates [L-b+1 .. L-a+1]; the alignment's 5' end is the high coordinate, so report
    // i = L-a+1 (start, high) and j = L-b+1 (stop, low), giving i > j as cmsearch does.
    let mut rc = dsq.to_vec();
    cm.abc.revcomp(&mut rc).expect("revcomp digital sequence");
    for h in scan_core::<S>(cm, &rc, w_max, Some(cutoff_bits)).1 {
        hits.push(Hit {
            score: h.score,
            evalue: evalue(cm, mode, h.score, searched),
            i: l - h.i + 1,
            j: l - h.j + 1,
            strand: Strand::Minus,
        });
    }

    // Best-first: by E-value ascending when available, otherwise by score descending.
    hits.sort_by(|a, b| match (a.evalue, b.evalue) {
        (Some(ea), Some(eb)) => ea
            .partial_cmp(&eb)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.score.total_cmp(&a.score)),
        _ => b.score.total_cmp(&a.score),
    });
    hits
}

/// Single-best convenience wrapper over [`scan_core`] (no multi-hit gamma resolution).
fn scan<S: Semiring>(cm: &Cm, dsq: &[Dsq], w_max: usize) -> CykScan {
    scan_core::<S>(cm, dsq, w_max, None).0
}

/// Shared scanning recursion, parameterized by the within-cell ⊕ ([`Semiring`]). Handles
/// both glocal and local mode (the latter when the CM was `configure_local`-ized): local
/// adds EL local-end initialization and ROOT local begins.
///
/// Returns the single best hit (`CykScan`) and, when `gamma_cutoff` is `Some`, the full set
/// of non-overlapping above-cutoff hits resolved by the gamma semi-HMM (in this strand's
/// local coordinates). With `gamma_cutoff = None` the gamma pass is skipped and the hit
/// vector is empty.
fn scan_core<S: Semiring>(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    gamma_cutoff: Option<f32>,
) -> (CykScan, Vec<RawHit>) {
    let l = dsq.len().saturating_sub(2);
    let w = w_max.min(l);
    let mut bestsc_per_j = vec![IMPOSSIBLE; l + 1];
    if l == 0 {
        return (
            CykScan {
                best: None,
                bestsc_per_j,
            },
            Vec::new(),
        );
    }
    let mut gamma = gamma_cutoff.map(|c| Gamma::new(l, c));

    let oesc = Oesc::build(cm);
    let kp = oesc.kp;
    let m = cm.m;

    // Compact indexing for BEGL_S states (the only states needing the W+1 rolling deck).
    let mut begl_idx = vec![usize::MAX; m];
    let mut n_begl = 0;
    for (v, slot) in begl_idx.iter_mut().enumerate() {
        if cm.stid[v] == stid::BEGL_S {
            *slot = n_begl;
            n_begl += 1;
        }
    }

    let width = w + 1;
    let mut alpha = vec![vec![IMPOSSIBLE; m * width]; 2];
    let mut alpha_begl = vec![vec![IMPOSSIBLE; n_begl.max(1) * width]; width];

    // Local-mode initialization: a state may end locally (v -> EL), emitting `dd` further
    // residues at `el_selfsc` each. init_sc(v, dd) = endsc[v] + el_selfsc*dd (IMPOSSIBLE in
    // glocal, where endsc is IMPOSSIBLE). Mirrors FCalcInitDPScores.
    let el = cm.el_selfsc;
    let init_sc = |v: usize, dd: usize| -> f32 {
        let e = cm.endsc[v];
        if e <= IMPOSSIBLE {
            IMPOSSIBLE
        } else {
            e + el * dd as f32
        }
    };
    // Local begin target states (first state of each MATP/MATL/MATR/BIF node). Empty in glocal.
    let begin_states: Vec<usize> = if cm.is_local {
        (1..m).filter(|&y| cm.beginsc[y] > IMPOSSIBLE).collect()
    } else {
        Vec::new()
    };

    // --- d=0 base case (empty-fragment scores), per cm_scan_mx_InitializeFloats. ---
    // Only E/S/D/B states get a finite d=0; emitters stay IMPOSSIBLE. Glocal => endsc is
    // IMPOSSIBLE everywhere (no local ends).
    for v in (0..m).rev() {
        if cm.stid[v] == stid::BEGL_S {
            let ych = cm.cfirst[v] as usize;
            let mut s = init_sc(v, 0); // endsc[v] (IMPOSSIBLE in glocal)
            for yoff in 0..cm.cnum[v] as usize {
                s = S::or(s, add(alpha[0][(ych + yoff) * width], cm.tsc[v][yoff]));
            }
            let bi = begl_idx[v];
            for deck in &mut alpha_begl {
                deck[bi * width] = s;
            }
        } else if cm.sttype[v] == st::E {
            alpha[0][v * width] = 0.0;
            alpha[1][v * width] = 0.0;
        } else if cm.sttype[v] == st::S || cm.sttype[v] == st::D {
            let ych = cm.cfirst[v] as usize;
            let mut s = init_sc(v, 0);
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
                    let mut sc = init_sc(v, d); // B sd=0; IMPOSSIBLE unless v can local-end
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
                let cd = d - sd;
                let mut sc = init_sc(v, cd); // local-end path: emit sd, then EL emits cd
                for yoff in 0..cnum {
                    sc = S::or(
                        sc,
                        add(alpha[child_row][(ych + yoff) * width + cd], tsc[yoff]),
                    );
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

        // ROOT_S (v = 0). Glocal: transitions to children. Local: ROOT transitions are
        // zeroed and entry is via local begins into internal start states. Report hits d>=1.
        let tsc0 = &cm.tsc[0];
        let ych = cm.cfirst[0] as usize;
        let cnum0 = cm.cnum[0] as usize;
        let mut bestsc_j = IMPOSSIBLE;
        // Gamma cell init (UpdateGammaHitMx): default to "no hit ends at j".
        if let Some(g) = gamma.as_mut() {
            g.mx[j] = g.mx[j - 1];
            g.gback[j] = 0;
            g.savesc[j] = IMPOSSIBLE;
        }
        for d in 1..=dx {
            let mut sc = IMPOSSIBLE;
            for yoff in 0..cnum0 {
                sc = S::or(sc, add(alpha[cur][(ych + yoff) * width + d], tsc0[yoff]));
            }
            // Local begins: ROOT_S -> internal start state y (never a BEGL_S).
            for &y in &begin_states {
                sc = S::or(sc, add(alpha[cur][y * width + d], cm.beginsc[y]));
            }
            alpha[cur][d] = sc; // ROOT stored at v=0 offset (read by nobody)
            if sc > bestsc_j {
                bestsc_j = sc;
            }
            if best.is_none_or(|b| sc > b.score) {
                best = Some(CykHit {
                    score: sc,
                    i: j - d + 1,
                    j,
                });
            }
            // Gamma update: a hit of length d ending at j, chained to the best resolution
            // ending before its start. Keep it if it improves the cumulative total.
            if let Some(g) = gamma.as_mut() {
                if sc > IMPOSSIBLE {
                    let i = j - d + 1;
                    let cumulative = g.mx[i - 1] + sc;
                    if cumulative > g.mx[j] {
                        g.mx[j] = cumulative;
                        g.gback[j] = i;
                        g.savesc[j] = sc;
                    }
                }
            }
        }
        bestsc_per_j[j] = bestsc_j;
    }

    let hits = gamma.map(|g| g.traceback(l)).unwrap_or_default();
    (CykScan { best, bestsc_per_j }, hits)
}
