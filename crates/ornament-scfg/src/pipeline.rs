//! Accelerated CM search: a cheap p7 filter cascade selects windows, the CM scans only the
//! survivors. A faithful-in-spirit port of `cm_pipeline.c`'s structure (p7 filter sweep →
//! surviving windows → CM stages), with an MSV → Forward filter cascade.
//!
//! Why this is correct (no true hits dropped vs the brute-force scan):
//!
//! 1. **Containment.** The strand is tiled into windows of length `2W` stepping by `W`
//!    (overlap `W`). Every CM hit has length `≤ W` (the model's max hit width), so any hit
//!    is fully contained in at least one tile — a length-`≤W` interval always fits inside
//!    one of these overlapping `2W` tiles.
//! 2. **Filtering (cheapest-first cascade).** Each tile is first scored by the reduced-precision
//!    uint8 **MSV** prefilter (16-lane SIMD); only tiles passing MSV's `F1` pay for the striped
//!    **Forward** filter, kept iff their Forward P-value passes `F3` (HMMER's third threshold,
//!    `1e-5`). MSV's `F1` is set with a margin so it keeps every tile Forward would (a real
//!    homolog scores far above either threshold), so the cascade prunes only signal-free tiles.
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

use ornament_alphabet::Dsq;
use ornament_hmm::profile::{bg_freqs, P7Profile};
use ornament_hmm::{
    forward_pvalue, msv_bits, msv_pvalue, parse_p7_hmm, viterbi_bits, viterbi_pvalue, P7Hmm,
};
use rayon::prelude::*;

use crate::evalues::{evalue, SearchMode};
use crate::model::Cm;
use crate::qdb::{calc_qdb_bands, QdbBands};
use crate::search::{scan_subseq, Hit, Strand};
use crate::InfernalError;

/// Tunable thresholds for the windowing pipeline. Defaults are HMMER's full p7 filter cascade
/// (`cm_pipeline.c`): MSV `F1` → Viterbi `F2` → Forward `F3`, so the native pipeline does the
/// same filtering work as `cmsearch` out of the box.
#[derive(Debug, Clone, Copy)]
pub struct PipelineParams {
    /// MSV-filter P-value threshold (HMMER's `F1`). A tile survives to the Viterbi stage iff
    /// its MSV P-value `≤ f1`. Default `0.02` (HMMER's default).
    pub f1: f64,
    /// Viterbi-filter P-value threshold (HMMER's `F2`). An MSV survivor survives to the Forward
    /// stage iff its Viterbi P-value `≤ f2`. Default `1e-3` (HMMER's default).
    pub f2: f64,
    /// Forward-filter P-value threshold (HMMER's `F3`). A Viterbi survivor survives to the CM
    /// stage iff its Forward P-value `≤ f3`. Default `1e-5`.
    pub f3: f64,
    /// Final CM stage: Inside (sum-over-parses) when `true`, CYK (best-parse) when `false`.
    pub do_inside: bool,
    /// Offload the MSV prefilter to the GPU (batch all tiles of a strand in one kernel) when the
    /// crate is built with the `cuda` feature and a device is available. A no-op otherwise, and
    /// it never changes the reported hits: GPU MSV is bit-parity-equal to the CPU MSV filter, so
    /// the survivor set (and thus the CM stage) is identical — only the prefilter runs faster.
    pub gpu_msv: bool,
}

impl Default for PipelineParams {
    fn default() -> Self {
        PipelineParams {
            f1: MSV_F1,
            f2: VIT_F2,
            f3: 1e-5,
            do_inside: false,
            gpu_msv: false,
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
    if l == 0 {
        return Ok((Vec::new(), stats));
    }

    // Configure the p7 Forward filter once. On x86_64 this builds the striped-SSE odds-space
    // profile (the per-tile sweep clones it cheaply and resets only the length model); the
    // length-independent base is shared read-only across both strands.
    let hmm = parse_filter_hmm(cm)?;
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(&hmm, cm.abc.as_ref(), &bg, l);
    let nj = prof.nj;
    let filt = fwdfilter::build(&prof);
    let msv = msvfilter::build(&prof); // cheap uint8 MSV prefilter (built once)
    let vit = vitfilter::build(&prof); // striped Viterbi filter (middle cascade stage, built once)
                                       // Query-dependent bands for the CM stage, computed once and shared across every survivor
                                       // window (and both strands). A safe `beta` keeps every real hit's optimal parse in-band, so
                                       // the banded scan reports identical hits to the unbanded one — just faster.
    let bands = calc_qdb_bands(cm, QdbBands::DEFAULT_BETA);
    let w = w_max.min(l).max(1);
    let searched = 2.0 * l as f64;
    let do_inside = params.do_inside;

    // GPU MSV prefilter: on only when built with `cuda` and requested (via the param or the
    // `ORNAMENT_GPU=1` env override). Without the feature it's forced off; `strand_survivors`
    // additionally falls back to the CPU filter if no device is available at runtime.
    #[cfg(feature = "cuda")]
    let use_gpu = params.gpu_msv || std::env::var("ORNAMENT_GPU").as_deref() == Ok("1");
    #[cfg(not(feature = "cuda"))]
    let use_gpu = {
        let _ = params.gpu_msv;
        false
    };

    // Reverse-complement strand. A window survivor in revcomp coordinate `x` maps back to
    // original coordinate `L - x + 1`; the alignment's 5' end is the high coordinate, so a
    // hit `i..j` (i ≤ j in revcomp) is reported `i' > j'`, matching `cmsearch`.
    let mut rc = dsq.to_vec();
    cm.abc.revcomp(&mut rc)?;

    // Each strand filters its survivors then CM-scans each survivor window; the windows are
    // independent (a hit reads only its own span), so they're scanned in parallel and the two
    // strands run concurrently via `rayon::join`. Output is order-independent: hits are sorted
    // below, so the result is identical regardless of thread count.
    let plus = || {
        let (surv, nw) =
            strand_survivors(&msv, &vit, &filt, nj, &hmm, &prof, dsq, w, &params, use_gpu);
        let hits: Vec<Hit> = surv
            .par_iter()
            .flat_map_iter(|&(s, e)| {
                let sub = subseq(dsq, s, e);
                scan_subseq(cm, &sub, w_max, cutoff_bits, do_inside, Some(&bands))
                    .into_iter()
                    .map(move |(score, i, j)| Hit {
                        score,
                        evalue: evalue(cm, mode, score, searched),
                        i: i + s - 1,
                        j: j + s - 1,
                        strand: Strand::Plus,
                    })
            })
            .collect();
        (surv, nw, hits)
    };
    let minus = || {
        let (surv, nw) =
            strand_survivors(&msv, &vit, &filt, nj, &hmm, &prof, &rc, w, &params, use_gpu);
        let hits: Vec<Hit> = surv
            .par_iter()
            .flat_map_iter(|&(s, e)| {
                let sub = subseq(&rc, s, e);
                scan_subseq(cm, &sub, w_max, cutoff_bits, do_inside, Some(&bands))
                    .into_iter()
                    .map(move |(score, i, j)| Hit {
                        score,
                        evalue: evalue(cm, mode, score, searched),
                        i: l - (i + s - 1) + 1,
                        j: l - (j + s - 1) + 1,
                        strand: Strand::Minus,
                    })
            })
            .collect();
        (surv, nw, hits)
    };
    let ((fsurv, fnw, fhits), (rsurv, rnw, rhits)) = rayon::join(plus, minus);

    stats.n_windows = fnw + rnw;
    stats.n_survivors = fsurv.len() + rsurv.len();
    stats.residues_to_cm = fsurv
        .iter()
        .chain(rsurv.iter())
        .map(|&(s, e)| e - s + 1)
        .sum();

    let mut hits: Vec<Hit> = fhits;
    hits.extend(rhits);

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

/// HMMER's default p7 filter cascade thresholds (P-value), for the prune-rate measurement.
const MSV_F1: f64 = 0.02; // MSV pass threshold (`--F1`)
const VIT_F2: f64 = 1e-3; // Viterbi pass threshold (`--F2`)

/// Per-stage tile pass counts for the p7 filter cascade (both strands summed). A measurement
/// used to size the value of an MSV → Viterbi → Forward cascade: at HMMER's default
/// thresholds, how many of the pipeline's tiles would each filter let through? It does **not**
/// change search output — it is pure instrumentation over the same tiling
/// [`cm_pipeline_search`] uses.
#[derive(Debug, Clone, Copy, Default)]
pub struct FilterPruneStats {
    /// Total tiles scored (both strands).
    pub n_tiles: usize,
    /// Tiles passing MSV (`P_MSV ≤ 0.02`).
    pub n_msv_pass: usize,
    /// Tiles passing Viterbi (`P_Vit ≤ 1e-3`).
    pub n_vit_pass: usize,
    /// Tiles passing Forward (`P_Fwd ≤ params.f3`) — identical to the live pipeline's survivor
    /// tiles (pre-merge), since it reuses the same Forward decision.
    pub n_fwd_pass: usize,
}

/// Measure the per-stage prune rates of an MSV → Viterbi → Forward cascade over the pipeline's
/// tiling, to decide whether a cheaper-first cascade is worth building. The Forward column uses
/// the exact same striped-Forward decision as [`cm_pipeline_search`] (so `n_fwd_pass` is its
/// survivor-tile count); MSV/Viterbi use the scalar filters at HMMER's default thresholds.
///
/// This is a diagnostic, not part of search — it returns counts only and leaves output unchanged.
pub fn measure_filter_prune_rates(
    cm: &Cm,
    dsq: &[Dsq],
    w_max: usize,
    params: PipelineParams,
) -> Result<FilterPruneStats, InfernalError> {
    let l = dsq.len().saturating_sub(2);
    let mut stats = FilterPruneStats::default();
    if l == 0 {
        return Ok(stats);
    }

    let hmm = parse_filter_hmm(cm)?;
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(&hmm, cm.abc.as_ref(), &bg, l);
    let nj = prof.nj;
    let filt = fwdfilter::build(&prof);
    let abc = cm.abc.as_ref();
    let w = w_max.min(l).max(1);

    let mut rc = dsq.to_vec();
    cm.abc.revcomp(&mut rc)?;

    for seq in [dsq, rc.as_slice()] {
        let (nt, nm, nv, nf) = measure_strand(seq, w, &hmm, abc, &filt, nj, params.f3);
        stats.n_tiles += nt;
        stats.n_msv_pass += nm;
        stats.n_vit_pass += nv;
        stats.n_fwd_pass += nf;
    }
    Ok(stats)
}

/// Tile-by-tile MSV/Viterbi/Forward pass counts for one strand (returns
/// `(n_tiles, n_msv, n_vit, n_fwd)`). Tiling matches [`strand_survivors`] exactly.
#[allow(clippy::too_many_arguments)]
fn measure_strand(
    seq: &[Dsq],
    w: usize,
    hmm: &P7Hmm,
    abc: &ornament_alphabet::Alphabet,
    filt: &fwdfilter::Profile,
    nj: f32,
    f3: f64,
) -> (usize, usize, usize, usize) {
    let l = seq.len().saturating_sub(2);
    let tile = (2 * w).max(1);
    let step = w.max(1);
    let mut tiles: Vec<(usize, usize)> = Vec::new();
    let mut start = 1usize;
    loop {
        let end = (start + tile - 1).min(l);
        tiles.push((start, end));
        if end >= l {
            break;
        }
        start += step;
    }

    tiles
        .par_iter()
        .map(|&(s, e)| {
            let sub = subseq(seq, s, e);
            let msv = msv_pvalue(hmm, msv_bits(hmm, abc, &sub)) <= MSV_F1;
            let vit = viterbi_pvalue(hmm, viterbi_bits(hmm, abc, &sub)) <= VIT_F2;
            let fwd = forward_pvalue(hmm, fwdfilter::window_bits(filt, nj, seq, s, e)) <= f3;
            (1usize, msv as usize, vit as usize, fwd as usize)
        })
        .reduce(
            || (0, 0, 0, 0),
            |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3),
        )
}

/// Parse the CM's embedded HMMER3/f filter HMM (`fp7`), erroring if it is absent.
fn parse_filter_hmm(cm: &Cm) -> Result<P7Hmm, InfernalError> {
    let text = cm
        .fp7_text
        .as_deref()
        .ok_or_else(|| InfernalError::Parse("CM has no embedded p7 filter HMM (fp7)".into()))?;
    parse_p7_hmm(text).map_err(|e| InfernalError::Parse(format!("filter HMM: {e}")))
}

/// Filter one strand through the full MSV → Viterbi → Forward cascade and return the merged,
/// W-padded survivor windows (1-based inclusive `(start, end)` in this strand's coordinates)
/// plus the number of tiles scored.
///
/// Tiling is length `2W`, step `W` (overlap `W`) so every length-`≤W` hit is contained in
/// at least one tile. Survivors are padded by `W` and merged into maximal disjoint regions.
#[allow(clippy::too_many_arguments)]
fn strand_survivors(
    msv: &msvfilter::Profile,
    vit: &vitfilter::Profile,
    filt: &fwdfilter::Profile,
    nj: f32,
    hmm: &P7Hmm,
    prof: &P7Profile,
    dsq: &[Dsq],
    w: usize,
    params: &PipelineParams,
    use_gpu: bool,
) -> (Vec<(usize, usize)>, usize) {
    let l = dsq.len().saturating_sub(2);
    let tile = (2 * w).max(1);
    let step = w.max(1);

    // Enumerate the tile spans (1-based inclusive) left-to-right, deterministically.
    let mut tiles: Vec<(usize, usize)> = Vec::new();
    let mut start = 1usize;
    loop {
        let end = (start + tile - 1).min(l);
        tiles.push((start, end));
        if end >= l {
            break;
        }
        start += step;
    }

    // Optional GPU MSV prefilter: score every tile's MSV in one batched kernel and return the
    // per-tile pass mask (indexed by tile position). `None` ⇒ score MSV per-tile on the CPU below.
    // The mask is bit-parity-equal to the CPU filter, so survivors are unchanged; it just moves
    // the cheapest, highest-volume stage off the CPU. Falls back to CPU if no device is available.
    let msv_mask: Option<Vec<bool>> = gpu_msv_mask(prof, hmm, dsq, &tiles, params.f1, use_gpu);

    // Cascade each tile cheapest-first (HMMER's `cm_pipeline.c`): the uint8 MSV prunes the large
    // majority, its survivors pay for the striped Viterbi, and only *those* survivors pay for the
    // (costliest) striped Forward. Each task clones the length-independent profiles and resets
    // only their length model, so tiles run in parallel without contention; the order-preserving
    // collect keeps the subsequent merge (assumes non-decreasing `start`) deterministic.
    let (f1, f2, f3) = (params.f1, params.f2, params.f3);
    let kept: Vec<Option<(usize, usize)>> = tiles
        .par_iter()
        .enumerate()
        .map(|(ti, &(s, e))| {
            let msv_ok = match &msv_mask {
                Some(mask) => mask[ti],                             // GPU batch decision
                None => msvfilter::passes(msv, hmm, dsq, s, e, f1), // CPU per-tile
            };
            if !msv_ok {
                return None; // MSV prefilter pruned this tile
            }
            if viterbi_pvalue(hmm, vitfilter::window_bits(vit, nj, dsq, s, e)) > f2 {
                return None; // Viterbi filter pruned this tile
            }
            if forward_pvalue(hmm, fwdfilter::window_bits(filt, nj, dsq, s, e)) <= f3 {
                let ps = s.saturating_sub(w).max(1);
                let pe = (e + w).min(l);
                Some((ps, pe))
            } else {
                None
            }
        })
        .collect();

    let mut survivors: Vec<(usize, usize)> = Vec::new();
    for (ps, pe) in kept.into_iter().flatten() {
        push_merge(&mut survivors, ps, pe);
    }
    (survivors, tiles.len())
}

/// GPU MSV pass mask for `tiles` (1-based inclusive spans into `dsq`): one batched kernel over the
/// whole strand, returning per-tile `P_MSV ≤ f1`. `None` when GPU is off/unbuilt or unavailable at
/// runtime, so the caller uses the CPU per-tile filter. The mask matches the CPU MSV filter
/// bit-for-bit (same uint8 recurrence + `msv_pvalue`), so it never changes the survivor set.
#[cfg(feature = "cuda")]
fn gpu_msv_mask(
    prof: &P7Profile,
    hmm: &P7Hmm,
    dsq: &[Dsq],
    tiles: &[(usize, usize)],
    f1: f64,
    use_gpu: bool,
) -> Option<Vec<bool>> {
    use ornament_hmm::msv_pvalue;
    if !use_gpu || tiles.is_empty() {
        return None;
    }
    // Any failure (no device, alloc error) → None → graceful CPU fallback.
    let bp = ornament_gpu::ByteMsvProfile::new(prof);
    let strand = ornament_gpu::DeviceStrand::upload(dsq).ok()?;
    let ctx = ornament_gpu::GpuContext::new().ok()?;
    let gtiles = ornament_gpu::Tiles::from_spans(tiles);
    let raw = ctx.msv_u8_nats(&bp, &strand, &gtiles).ok()?;
    Some(
        tiles
            .iter()
            .zip(raw)
            .map(|(&(s, e), r)| {
                let bits = ornament_gpu::nats_to_bits(r, e - s + 1);
                msv_pvalue(hmm, bits) <= f1
            })
            .collect(),
    )
}

/// No-CUDA stub: always `None`, so the pipeline uses the CPU MSV filter.
#[cfg(not(feature = "cuda"))]
fn gpu_msv_mask(
    _prof: &P7Profile,
    _hmm: &P7Hmm,
    _dsq: &[Dsq],
    _tiles: &[(usize, usize)],
    _f1: f64,
    _use_gpu: bool,
) -> Option<Vec<bool>> {
    None
}

/// The p7 Forward filter, abstracted over its implementation: the striped-SSE odds-space
/// kernel on x86_64, the scalar log-space DP elsewhere. `build` makes the length-independent
/// profile once; `window_bits` scores `dsq[start..=end]` (cloning the profile and resetting
/// only its length model so tiles can run in parallel) and returns the Forward **bit** score,
/// matching [`ornament_hmm::forward_bits`].
#[cfg(target_arch = "x86_64")]
mod fwdfilter {
    use super::{subseq, Dsq, P7Profile};
    use ornament_hmm::profile::null_one;

    pub(super) type Profile = ornament_hmm::StripedProfile;

    pub(super) fn build(prof: &P7Profile) -> Profile {
        ornament_hmm::StripedProfile::new(prof)
    }

    pub(super) fn window_bits(
        base: &Profile,
        nj: f32,
        dsq: &[Dsq],
        start: usize,
        end: usize,
    ) -> f32 {
        let len = end - start + 1;
        let sub = subseq(dsq, start, end);
        let mut sp = base.clone();
        sp.set_length(len, nj);
        (sp.score_nats(&sub) - null_one(len)) / std::f32::consts::LN_2
    }
}

#[cfg(not(target_arch = "x86_64"))]
mod fwdfilter {
    use super::{subseq, Dsq, P7Profile};
    use ornament_hmm::forward_nats;
    use ornament_hmm::profile::null_one;

    pub(super) type Profile = P7Profile;

    pub(super) fn build(prof: &P7Profile) -> Profile {
        prof.clone()
    }

    pub(super) fn window_bits(
        base: &Profile,
        _nj: f32,
        dsq: &[Dsq],
        start: usize,
        end: usize,
    ) -> f32 {
        let len = end - start + 1;
        let sub = subseq(dsq, start, end);
        let mut p = base.clone();
        p.reconfig_length(len);
        (forward_nats(&p, &sub) - null_one(len)) / std::f32::consts::LN_2
    }
}

/// The p7 Viterbi filter (middle cascade stage), abstracted over its implementation: the
/// striped-SIMD max-plus kernel on x86_64, the scalar log-space DP elsewhere. `build` makes the
/// length-independent profile once; `window_bits` scores `dsq[start..=end]` (cloning + resetting
/// only the length model so tiles run in parallel) and returns the Viterbi **bit** score,
/// matching [`ornament_hmm::viterbi_bits`].
#[cfg(target_arch = "x86_64")]
mod vitfilter {
    use super::{subseq, Dsq, P7Profile};
    use ornament_hmm::profile::null_one;

    pub(super) type Profile = ornament_hmm::StripedVitProfile;

    pub(super) fn build(prof: &P7Profile) -> Profile {
        ornament_hmm::StripedVitProfile::new(prof)
    }

    pub(super) fn window_bits(
        base: &Profile,
        nj: f32,
        dsq: &[Dsq],
        start: usize,
        end: usize,
    ) -> f32 {
        let len = end - start + 1;
        let sub = subseq(dsq, start, end);
        let mut sp = base.clone();
        sp.set_length(len, nj);
        (sp.score_nats(&sub) - null_one(len)) / std::f32::consts::LN_2
    }
}

#[cfg(not(target_arch = "x86_64"))]
mod vitfilter {
    use super::{subseq, Dsq, P7Profile};
    use ornament_hmm::profile::null_one;
    use ornament_hmm::viterbi_nats;

    pub(super) type Profile = P7Profile;

    pub(super) fn build(prof: &P7Profile) -> Profile {
        prof.clone()
    }

    pub(super) fn window_bits(
        base: &Profile,
        _nj: f32,
        dsq: &[Dsq],
        start: usize,
        end: usize,
    ) -> f32 {
        let len = end - start + 1;
        let sub = subseq(dsq, start, end);
        let mut p = base.clone();
        p.reconfig_length(len);
        (viterbi_nats(&p, &sub) - null_one(len)) / std::f32::consts::LN_2
    }
}

/// The reduced-precision (uint8) MSV prefilter — the cheapest cascade stage. On x86_64 it's the
/// 16-lane striped int8 MSV; elsewhere it's a no-op pass-through (the Forward filter still gates).
/// `passes` is whether a tile clears the MSV P-value threshold `f1`; a tile that fails is pruned
/// without ever running the costlier Forward DP. Strong hits saturate the uint8 range → `+inf`
/// score → always pass, so the prefilter cannot drop a real hit (verified by pipeline parity).
#[cfg(target_arch = "x86_64")]
mod msvfilter {
    use super::{subseq, Dsq, P7Hmm, P7Profile};
    use ornament_hmm::{msv_pvalue, MsvProfile};

    pub(super) type Profile = MsvProfile;

    pub(super) fn build(prof: &P7Profile) -> Profile {
        MsvProfile::new(prof)
    }

    pub(super) fn passes(
        base: &Profile,
        hmm: &P7Hmm,
        dsq: &[Dsq],
        start: usize,
        end: usize,
        f1: f64,
    ) -> bool {
        let mut sp = base.clone();
        sp.set_length(end - start + 1);
        let sub = subseq(dsq, start, end);
        msv_pvalue(hmm, sp.score_bits(&sub)) <= f1
    }
}

#[cfg(not(target_arch = "x86_64"))]
mod msvfilter {
    use super::{Dsq, P7Hmm, P7Profile};

    pub(super) type Profile = ();

    pub(super) fn build(_prof: &P7Profile) -> Profile {}

    pub(super) fn passes(
        _base: &Profile,
        _hmm: &P7Hmm,
        _dsq: &[Dsq],
        _start: usize,
        _end: usize,
        _f1: f64,
    ) -> bool {
        true // no MSV prefilter off x86; the Forward filter still gates every tile
    }
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
