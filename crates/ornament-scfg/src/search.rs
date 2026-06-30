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

use ornament_alphabet::Dsq;

use crate::config::IMPOSSIBLE;
use crate::emit::Oesc;
use crate::evalues::{evalue, SearchMode};
use crate::model::{emits_right, n_emit, st, stid, Cm};
use crate::qdb::QdbBands;

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

    /// Fold a child transition into the accumulator row: `acc[i] ⊕= (src[i] ⊗ tsc)` for every
    /// `i in 0..acc.len()` (`src` must be at least that long). This is the hot inner step of
    /// the scanning DP — the `cnum`-way child max/sum over a contiguous `d`-slice. The default
    /// is the scalar fold via [`Semiring::or`] and the sentinel-clamping [`add`]; [`MaxPlus`]
    /// overrides it with a branchless max-plus loop the compiler auto-vectorizes (AVX2
    /// `vaddps`/`vmaxps`).
    #[inline]
    fn accumulate(acc: &mut [f32], src: &[f32], tsc: f32) {
        for (a, &s) in acc.iter_mut().zip(src) {
            *a = Self::or(*a, add(s, tsc));
        }
    }
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

    /// Branchless max-plus fold. Bit-identical to the default scalar path despite skipping the
    /// sentinel-clamping `add`: `src[i] + tsc` adds at most **one** `IMPOSSIBLE` (the child may
    /// be impossible; `tsc` is a finite transition score, or itself `IMPOSSIBLE`), and
    /// `IMPOSSIBLE` (`-1e36`) so dwarfs a finite addend that the f32 sum stays exactly
    /// `IMPOSSIBLE`. The running `max` never drops below the caller's `>= IMPOSSIBLE` init, so a
    /// two-sentinel `-2e36` can never become the chosen value. No `or`/`add` calls and no
    /// sentinel branch ⇒ the loop vectorizes.
    #[inline]
    fn accumulate(acc: &mut [f32], src: &[f32], tsc: f32) {
        // Runtime dispatch: AVX2 256-bit `vmaxps` when available (the whole cluster has it via
        // the x86-64-v3 baseline), else the scalar fallback. Both are bit-identical — AVX2
        // reserved here for the explicit kernel; AVX-512 could slot in the same way.
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                // SAFETY: guarded by the feature detection above.
                unsafe { accumulate_max_avx2(acc, src, tsc) };
                return;
            }
        }
        accumulate_max_scalar(acc, src, tsc);
    }
}

/// Scalar max-plus child fold: `acc[i] = max(acc[i], src[i] + tsc)`. The fallback for
/// [`accumulate_max_avx2`] and the reference the SIMD path must match bit-for-bit.
#[inline]
fn accumulate_max_scalar(acc: &mut [f32], src: &[f32], tsc: f32) {
    let n = acc.len();
    let src = &src[..n];
    for i in 0..n {
        // `a < v ? v : a` mirrors `vmaxps(a, v)` exactly (no NaN/-0 here).
        let a = acc[i];
        let v = src[i] + tsc;
        acc[i] = if a < v { v } else { a };
    }
}

/// AVX2 max-plus child fold (8 lanes via `vmaxps`). Bit-identical to [`accumulate_max_scalar`]:
/// `_mm256_max_ps(a, v)` returns `a < v ? v : a` lane-wise, and the scalar tail uses the same
/// rule. `src` is guaranteed at least `acc.len()` long by the caller.
///
/// # Safety
/// Requires the `avx2` target feature (checked by the caller via `is_x86_feature_detected!`).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn accumulate_max_avx2(acc: &mut [f32], src: &[f32], tsc: f32) {
    use std::arch::x86_64::*;
    let n = acc.len();
    let tv = _mm256_set1_ps(tsc);
    let (ap, sp) = (acc.as_mut_ptr(), src.as_ptr());
    let mut i = 0;
    while i + 8 <= n {
        let a = _mm256_loadu_ps(ap.add(i));
        let v = _mm256_add_ps(_mm256_loadu_ps(sp.add(i)), tv);
        _mm256_storeu_ps(ap.add(i), _mm256_max_ps(a, v));
        i += 8;
    }
    while i < n {
        let a = *ap.add(i);
        let v = *sp.add(i) + tsc;
        *ap.add(i) = if a < v { v } else { a };
        i += 1;
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
/// Unbanded — the brute-force reference; for the faster banded form use [`cyk_search_banded`].
pub fn cyk_search(cm: &Cm, dsq: &[Dsq], w_max: usize, cutoff_bits: f32) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalCyk
    } else {
        SearchMode::GlocalCyk
    };
    search::<MaxPlus>(cm, dsq, w_max, cutoff_bits, mode, None)
}

/// Query-dependent-banded [`cyk_search`]. Identical hits to the unbanded form (the bands keep
/// every real hit's optimal parse in range) at a fraction of the DP cost. Compute `bands` once
/// per model via [`crate::qdb::calc_qdb_bands`] and reuse across calls.
pub fn cyk_search_banded(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    bands: &QdbBands,
) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalCyk
    } else {
        SearchMode::GlocalCyk
    };
    search::<MaxPlus>(cm, dsq, w_max, cutoff_bits, mode, Some(bands))
}

/// Inside counterpart of [`cyk_search`] (sum-over-parses scoring). Unbanded.
pub fn inside_search(cm: &Cm, dsq: &[Dsq], w_max: usize, cutoff_bits: f32) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalInside
    } else {
        SearchMode::GlocalInside
    };
    search::<LogSum>(cm, dsq, w_max, cutoff_bits, mode, None)
}

/// Query-dependent-banded [`inside_search`].
pub fn inside_search_banded(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    bands: &QdbBands,
) -> Vec<Hit> {
    let mode = if cm.is_local {
        SearchMode::LocalInside
    } else {
        SearchMode::GlocalInside
    };
    search::<LogSum>(cm, dsq, w_max, cutoff_bits, mode, Some(bands))
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
    bands: Option<&QdbBands>,
) -> Vec<Hit> {
    let l = dsq.len().saturating_sub(2);
    let searched = 2.0 * l as f64;

    // The two strands are independent DP scans; run them concurrently. `rayon::join` is a
    // synchronous fork-join, so `rc` (borrowed by the second closure) outlives the call.
    let mut rc = dsq.to_vec();
    cm.abc.revcomp(&mut rc).expect("revcomp digital sequence");
    let (fwd, rev) = rayon::join(
        || scan_core::<S>(cm, dsq, w_max, Some(cutoff_bits), bands).1,
        || scan_core::<S>(cm, &rc, w_max, Some(cutoff_bits), bands).1,
    );

    let mut hits: Vec<Hit> = Vec::with_capacity(fwd.len() + rev.len());

    // Forward strand: coordinates are already original.
    for h in fwd {
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
    for h in rev {
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
    scan_core::<S>(cm, dsq, w_max, None, None).0
}

/// Multi-hit gamma scan over a *standalone* sub-sequence `sub` (its own sentinels at 0 and
/// `len+1`), returning each above-`cutoff_bits` hit as `(score, i, j)` in `sub`'s local
/// 1-based coordinates. `do_inside` selects Inside (sum-over-parses) vs CYK (best-parse).
///
/// This is the primitive the [`crate::pipeline`] windowing path runs on filter survivors.
/// It reuses [`scan_core`] unchanged, so scores and coordinates are bit-identical to the
/// brute-force whole-strand scan: a CM hit rooted at ROOT_S over `i..j` reads only
/// `sub[i..=j]`, so flanking context outside a window can never change its parse — provided
/// the window fully contains the hit (the caller guarantees this with W-padding).
pub(crate) fn scan_subseq(
    cm: &Cm,
    sub: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    do_inside: bool,
    bands: Option<&QdbBands>,
) -> Vec<(f32, usize, usize)> {
    let raw = if do_inside {
        scan_core::<LogSum>(cm, sub, w_max, Some(cutoff_bits), bands).1
    } else {
        scan_core::<MaxPlus>(cm, sub, w_max, Some(cutoff_bits), bands).1
    };
    raw.into_iter().map(|h| (h.score, h.i, h.j)).collect()
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
    bands: Option<&QdbBands>,
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
    // Reused per-state accumulator for the vectorized child-transition fold (holds the running
    // ⊕ over children for each `d`, before the per-`d` emission is applied).
    let mut acc = vec![IMPOSSIBLE; width];

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

            // QDB band for this state: only `d in [dmin_v, dmax_v]` is computed/stored. The band
            // is j-independent, so out-of-band cells are never written and stay at their
            // one-time IMPOSSIBLE init — reads of them are correctly IMPOSSIBLE. Unbanded =
            // the full `[0, dx]` range (identical to the pre-QDB scan).
            let (dmin_v, dmax_v) = bands.map_or((0, dx), |b| (b.dmin[v], b.dmax[v].min(dx)));

            if sttype == st::B {
                // Bifurcation: sc[d] = ⊕_{k=0..d} left(j-k, d-k) ⊗ right(j, k). Reordering to
                // outer-`k` makes the left child's cells contiguous in `d-k`, so each split is
                // the same vectorized child-fold as the emitting states (`S::accumulate`), and
                // the rolling-deck modulo is computed once per `k` instead of per `(d, k)`.
                let wbi = begl_idx[cm.cfirst[v] as usize];
                let ych = cm.cnum[v] as usize;
                let base = wbi * width;
                let dlo = dmin_v.max(1);
                let dhi = dmax_v;
                if dhi < dlo {
                    continue;
                }
                // Children's bands: `k` = right-fragment length (BEGR via cnum), `d-k` = left
                // (BEGL via cfirst). Restricting both keeps the convolution in-band (Infernal's
                // kmin/kmax) — out-of-band child cells are IMPOSSIBLE anyway, so this is purely
                // skipping guaranteed-no-op work.
                let (lmin, lmax, rmin, rmax) = bands.map_or((0, dx, 0, dx), |b| {
                    let l = cm.cfirst[v] as usize;
                    (
                        b.dmin[l],
                        b.dmax[l].min(dx),
                        b.dmin[ych],
                        b.dmax[ych].min(dx),
                    )
                });
                #[allow(clippy::needless_range_loop)] // d feeds init_sc and indexes acc
                for d in dlo..=dhi {
                    acc[d] = init_sc(v, d); // B sd=0; IMPOSSIBLE unless v can local-end
                }
                for k in rmin..=rmax {
                    let right_k = alpha[cur][ych * width + k];
                    if right_k <= IMPOSSIBLE {
                        continue; // left ⊗ IMPOSSIBLE contributes nothing to any d
                    }
                    let deck = (j - k) % width;
                    // d in band, d >= k, and left length d-k in [lmin, lmax].
                    let d_start = dlo.max(k + lmin);
                    let d_end = dhi.min(k + lmax);
                    if d_end < d_start {
                        continue;
                    }
                    let n = d_end - d_start + 1;
                    let loff = base + (d_start - k); // left index (d-k) at d = d_start
                    S::accumulate(
                        &mut acc[d_start..=d_end],
                        &alpha_begl[deck][loff..loff + n],
                        right_k,
                    );
                }
                #[allow(clippy::needless_range_loop)] // d indexes acc and the alpha row offset
                for d in dlo..=dhi {
                    alpha[cur][v * width + d] = acc[d];
                }
                continue;
            }

            let ych = cm.cfirst[v] as usize;
            let cnum = cm.cnum[v] as usize;
            let tsc = &cm.tsc[v];
            let child_row = if right { prv } else { cur };
            let single = oesc.single[v].as_deref();
            let pairtab = oesc.pair[v].as_deref();

            // Banded d-range for this emitter/D/S state: `d in [sd.max(dmin_v), dmax_v]`.
            let dlo = sd.max(dmin_v);
            let dhi = dmax_v;
            if dhi < dlo {
                continue; // empty band / too few residues; all cells stay IMPOSSIBLE
            }
            // A same-row self-transition (`v` among its own current-row children) makes the
            // `d`-recurrence loop-carried: `alpha[cur][v][d]` reads `alpha[cur][v][d-1]` computed
            // earlier this loop. That's the insert-left (IL) case — it can't be split into a
            // bulk child fold, so keep the sequential recurrence. IR self-loops too but reads the
            // *previous* row (`prv`), so it stays in the vectorized path below.
            let self_loop = child_row == cur && (ych..ych + cnum).contains(&v);

            // `d` starts at `sd.max(1)`: emitters need `d >= sd`, and the `d=0` cell of S/D
            // states is a constant set in the base case that must not be overwritten here.
            if self_loop {
                // IL self-loop reads `alpha[cur][v][d-1]` (its own earlier cell). At `d = dlo`
                // that read is out-of-band → IMPOSSIBLE, exactly as the minimal-length parse
                // requires; processing `d` ascending keeps the recurrence correct.
                for d in dlo.max(1)..=dhi {
                    let cd = d - sd;
                    let mut sc = init_sc(v, cd);
                    #[allow(clippy::needless_range_loop)] // yoff indexes both the child row and tsc
                    for yoff in 0..cnum {
                        sc = S::or(
                            sc,
                            add(alpha[child_row][(ych + yoff) * width + cd], tsc[yoff]),
                        );
                    }
                    // IL is the only self-loop here (left emitter); store with its left emission.
                    let val = add(sc, single.unwrap()[dsq[j - d + 1] as usize]);
                    alpha[cur][v * width + d] = val;
                }
                continue;
            }

            // Phase 1: child-transition fold over the banded `cd = d - sd` range (the hot,
            // vectorized inner step). `cd` runs cd0..=cd1, i.e. `d` runs dlo..=dhi.
            let cd0 = dlo - sd;
            let cd1 = dhi - sd;
            #[allow(clippy::needless_range_loop)] // cd feeds init_sc and indexes acc
            for cd in cd0..=cd1 {
                acc[cd] = init_sc(v, cd); // local-end path: emit sd, then EL emits cd (IMPOSSIBLE glocal)
            }
            #[allow(clippy::needless_range_loop)] // yoff indexes both the child row and tsc
            for yoff in 0..cnum {
                let base = (ych + yoff) * width;
                S::accumulate(
                    &mut acc[cd0..=cd1],
                    &alpha[child_row][base + cd0..base + cd1 + 1],
                    tsc[yoff],
                );
            }
            // Phase 2: add the per-`d` emission score and store (begl states have no emission).
            for d in dlo.max(1)..=dhi {
                let sc = acc[d - sd];
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

#[cfg(test)]
mod simd_tests {
    use super::*;

    /// The vectorized [`MaxPlus::accumulate`] must be **bit-for-bit** identical to the generic
    /// scalar fold (`or(acc, add(src, tsc))`) — the SIMD-vs-scalar parity invariant from #5.
    /// Inputs respect the DP invariant that stored cells are always `>= IMPOSSIBLE` (either a
    /// finite score or exactly the sentinel); `tsc` may additionally be `IMPOSSIBLE`.
    #[test]
    fn maxplus_accumulate_is_bit_exact_vs_scalar() {
        let mut s: u64 = 0xDEAD_BEEF_1234_5678;
        let mut next = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s
        };
        // A value respecting `>= IMPOSSIBLE`: ~1/4 of the time the sentinel, else a finite
        // log-odds-scale score in roughly [-80, 40).
        let cell = |r: u64| -> f32 {
            if r & 3 == 0 {
                IMPOSSIBLE
            } else {
                ((r >> 8) % 12000) as f32 / 100.0 - 80.0
            }
        };

        for _ in 0..2000 {
            let n = 1 + (next() as usize % 40); // odd lengths exercise the SIMD remainder
            let cnum = 1 + (next() as usize % 4);
            let init: Vec<f32> = (0..n).map(|_| cell(next())).collect();
            let children: Vec<Vec<f32>> = (0..cnum)
                .map(|_| (0..n).map(|_| cell(next())).collect())
                .collect();
            // tsc: finite, or IMPOSSIBLE for an impossible transition.
            let tsc: Vec<f32> = (0..cnum)
                .map(|_| {
                    let r = next();
                    if r & 7 == 0 {
                        IMPOSSIBLE
                    } else {
                        (r >> 8) as i32 as f32 / 1e7
                    }
                })
                .collect();

            let mut got = init.clone();
            let mut want = init.clone();
            for y in 0..cnum {
                MaxPlus::accumulate(&mut got, &children[y], tsc[y]);
                for i in 0..n {
                    want[i] = MaxPlus::or(want[i], add(children[y][i], tsc[y]));
                }
            }
            for i in 0..n {
                assert_eq!(
                    got[i].to_bits(),
                    want[i].to_bits(),
                    "mismatch at i={i}: vectorized {} vs scalar {}",
                    got[i],
                    want[i]
                );
            }
        }
    }
}
