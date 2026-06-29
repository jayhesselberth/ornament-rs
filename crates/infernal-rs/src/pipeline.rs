//! Accelerated CM search: the p7 Forward filter selects windows, the CM scans only the
//! survivors. A faithful-in-spirit port of `cm_pipeline.c`'s structure (p7 filter sweep →
//! surviving windows → CM stages), reduced to a single Forward filter stage.
//!
//! Why this is correct (no true hits dropped vs the brute-force scan):
//!
//! 1. **Containment.** The strand is tiled into windows of length `2W` stepping by `W`
//!    (overlap `W`). Every CM hit has length `≤ W` (the model's max hit width), so any hit
//!    is fully contained in at least one tile — a length-`≤W` interval always fits inside
//!    one of these overlapping `2W` tiles.
//! 2. **Filtering.** Each tile is scored with p7 Forward and kept iff its Forward P-value
//!    passes `F3` (HMMER's third filter threshold, default `1e-5`). A real homolog scores
//!    tens of bits over its tile → P-value far below `1e-5`, so its containing tile always
//!    survives; only signal-free tiles are pruned.
//! 3. **Window-local CM = whole-strand CM.** A surviving tile is padded by `W` (so adjacent
//!    survivors merge and every contained hit keeps full flank room) and handed to the CM
//!    scan. Because a CM hit rooted at ROOT_S over `i..j` reads only `dsq[i..=j]`, scanning
//!    a window that *contains* the hit reproduces its score and coordinates exactly. Hits
//!    can't span the gaps between disjoint survivor regions (each fits in one tile ⊂ one
//!    region), so per-region gamma resolution equals the global one over the kept regions.
//!
//! The net effect: the CM DP — the dominant cost — runs only on survivor residues, while
//! the cheap O(W·M)-per-tile Forward sweep prunes the rest. Reported scores/E-values come
//! from the CM core exactly as in [`crate::search::cyk_search`].

use easel_rs::Dsq;
use hmmer_rs::profile::{bg_freqs, null_one, P7Profile};
use hmmer_rs::{forward_nats, forward_pvalue, parse_p7_hmm, P7Hmm};

use crate::evalues::{evalue, SearchMode};
use crate::model::Cm;
use crate::search::{scan_subseq, Hit, Strand};
use crate::InfernalError;

/// Tunable thresholds for the windowing pipeline.
#[derive(Debug, Clone, Copy)]
pub struct PipelineParams {
    /// Forward-filter P-value threshold (HMMER's `F3`). A window survives to the CM stage
    /// iff its Forward P-value `≤ f3`. Default `1e-5`.
    pub f3: f64,
    /// Final CM stage: Inside (sum-over-parses) when `true`, CYK (best-parse) when `false`.
    pub do_inside: bool,
}

impl Default for PipelineParams {
    fn default() -> Self {
        PipelineParams {
            f3: 1e-5,
            do_inside: false,
        }
    }
}

/// Bookkeeping on how much work the filter saved, summed over both strands.
#[derive(Debug, Clone, Copy, Default)]
pub struct PipelineStats {
    /// Total residues considered (both strands) = `2 * L`.
    pub residues_searched: usize,
    /// Residues actually handed to the CM DP (sum of survivor-window lengths, both strands).
    pub residues_to_cm: usize,
    /// Number of Forward-filter windows evaluated.
    pub n_windows: usize,
    /// Number of surviving (merged) CM windows.
    pub n_survivors: usize,
}

impl PipelineStats {
    /// Fraction of residues the filter forwarded to the CM (`1.0` = no pruning).
    pub fn cm_fraction(&self) -> f64 {
        if self.residues_searched == 0 {
            0.0
        } else {
            self.residues_to_cm as f64 / self.residues_searched as f64
        }
    }
}

/// Accelerated CYK search via the Forward-filter windowing pipeline. Drop-in for
/// [`crate::search::cyk_search`] (set `params.do_inside = false`), returning the same hits
/// above `cutoff_bits` plus a [`PipelineStats`] summary. Honours the CM's glocal/local mode.
pub fn cm_pipeline_search(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    params: PipelineParams,
) -> Result<(Vec<Hit>, PipelineStats), InfernalError> {
    let l = dsq.len().saturating_sub(2);
    let mode = match (cm.is_local, params.do_inside) {
        (false, false) => SearchMode::GlocalCyk,
        (false, true) => SearchMode::GlocalInside,
        (true, false) => SearchMode::LocalCyk,
        (true, true) => SearchMode::LocalInside,
    };

    let mut stats = PipelineStats {
        residues_searched: 2 * l,
        ..Default::default()
    };
    let mut hits: Vec<Hit> = Vec::new();
    if l == 0 {
        return Ok((hits, stats));
    }

    // Configure the p7 Forward filter once; the length model is reconfigured per window.
    let hmm = parse_filter_hmm(cm)?;
    let bg = bg_freqs(hmm.k);
    let mut prof = P7Profile::config_local(&hmm, cm.abc.as_ref(), &bg, l);
    let w = w_max.min(l).max(1);
    let searched = 2.0 * l as f64;

    // Forward strand: filter survivors, then CM scan each (coords are already original).
    let (surv, nw) = strand_survivors(&mut prof, &hmm, dsq, w, params.f3);
    stats.n_windows += nw;
    for &(s, e) in &surv {
        stats.n_survivors += 1;
        stats.residues_to_cm += e - s + 1;
        let sub = subseq(dsq, s, e);
        for (score, i, j) in scan_subseq(cm, &sub, w_max, cutoff_bits, params.do_inside) {
            hits.push(Hit {
                score,
                evalue: evalue(cm, mode, score, searched),
                i: i + s - 1,
                j: j + s - 1,
                strand: Strand::Plus,
            });
        }
    }

    // Reverse-complement strand. A window survivor in revcomp coordinate `x` maps back to
    // original coordinate `L - x + 1`; the alignment's 5' end is the high coordinate, so a
    // hit `i..j` (i ≤ j in revcomp) is reported `i' > j'`, matching `cmsearch`.
    let mut rc = dsq.to_vec();
    cm.abc.revcomp(&mut rc)?;
    let (surv, nw) = strand_survivors(&mut prof, &hmm, &rc, w, params.f3);
    stats.n_windows += nw;
    for &(s, e) in &surv {
        stats.n_survivors += 1;
        stats.residues_to_cm += e - s + 1;
        let sub = subseq(&rc, s, e);
        for (score, i, j) in scan_subseq(cm, &sub, w_max, cutoff_bits, params.do_inside) {
            hits.push(Hit {
                score,
                evalue: evalue(cm, mode, score, searched),
                i: l - (i + s - 1) + 1,
                j: l - (j + s - 1) + 1,
                strand: Strand::Minus,
            });
        }
    }

    // Best-first, identical ordering to `search::search` (E-value asc, else score desc).
    hits.sort_by(|a, b| match (a.evalue, b.evalue) {
        (Some(ea), Some(eb)) => ea
            .partial_cmp(&eb)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.score.total_cmp(&a.score)),
        _ => b.score.total_cmp(&a.score),
    });
    Ok((hits, stats))
}

/// Inside counterpart of [`cm_pipeline_search`] (sum-over-parses final stage).
pub fn cm_pipeline_inside(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    cutoff_bits: f32,
    mut params: PipelineParams,
) -> Result<(Vec<Hit>, PipelineStats), InfernalError> {
    params.do_inside = true;
    cm_pipeline_search(cm, dsq, w_max, cutoff_bits, params)
}

/// Parse the CM's embedded HMMER3/f filter HMM (`fp7`), erroring if it is absent.
fn parse_filter_hmm(cm: &Cm) -> Result<P7Hmm, InfernalError> {
    let text = cm
        .fp7_text
        .as_deref()
        .ok_or_else(|| InfernalError::Parse("CM has no embedded p7 filter HMM (fp7)".into()))?;
    parse_p7_hmm(text).map_err(|e| InfernalError::Parse(format!("filter HMM: {e}")))
}

/// Forward-filter one strand and return the merged, W-padded survivor windows (1-based
/// inclusive `(start, end)` in this strand's coordinates) plus the number of tiles scored.
///
/// Tiling is length `2W`, step `W` (overlap `W`) so every length-`≤W` hit is contained in
/// at least one tile. Survivors are padded by `W` and merged into maximal disjoint regions.
fn strand_survivors(
    prof: &mut P7Profile,
    hmm: &P7Hmm,
    dsq: &[Dsq],
    w: usize,
    f3: f64,
) -> (Vec<(usize, usize)>, usize) {
    let l = dsq.len().saturating_sub(2);
    let tile = (2 * w).max(1);
    let step = w.max(1);
    let mut survivors: Vec<(usize, usize)> = Vec::new();
    let mut n_windows = 0usize;

    let mut start = 1usize;
    loop {
        let end = (start + tile - 1).min(l);
        n_windows += 1;
        if forward_pvalue(hmm, window_forward_bits(prof, dsq, start, end)) <= f3 {
            let ps = start.saturating_sub(w).max(1);
            let pe = (end + w).min(l);
            push_merge(&mut survivors, ps, pe);
        }
        if end >= l {
            break;
        }
        start += step;
    }
    (survivors, n_windows)
}

/// Forward **bit** score of `dsq[start..=end]` against the length-reconfigured profile,
/// matching [`hmmer_rs::forward_bits`] but reusing the prebuilt profile.
fn window_forward_bits(prof: &mut P7Profile, dsq: &[Dsq], start: usize, end: usize) -> f32 {
    let len = end - start + 1;
    let sub = subseq(dsq, start, end);
    prof.reconfig_length(len);
    (forward_nats(prof, &sub) - null_one(len)) / std::f32::consts::LN_2
}

/// Build a standalone digital sub-sequence for `dsq[start..=end]` with fresh sentinels at 0
/// and `len+1` (copied from the parent's leading sentinel), as the DP routines expect.
fn subseq(dsq: &[Dsq], start: usize, end: usize) -> Vec<Dsq> {
    let sentinel = dsq[0];
    let mut out = Vec::with_capacity(end - start + 3);
    out.push(sentinel);
    out.extend_from_slice(&dsq[start..=end]);
    out.push(sentinel);
    out
}

/// Insert `(s, e)` into the ascending, disjoint survivor list, coalescing any overlapping or
/// touching regions so each CM scan covers a maximal contiguous span exactly once.
fn push_merge(regions: &mut Vec<(usize, usize)>, s: usize, e: usize) {
    if let Some(last) = regions.last_mut() {
        // Tiles are visited left-to-right, so `s` is non-decreasing: merge with the tail.
        if s <= last.1 + 1 {
            last.1 = last.1.max(e);
            return;
        }
    }
    regions.push((s, e));
}
