//! Native CM search path (pure-Rust `ornament-scfg` engine).
//!
//! This is the default `scan` engine: a drop-in replacement for the `cmsearch`
//! subprocess that needs no external binary or C toolchain. It parses a `.cm`,
//! configures it in local mode (the `cmsearch` default), runs the native CYK
//! scanner over every FASTA record (both strands, multi-hit, gamma overlap
//! resolution), and yields the same [`CMHit`] records the CLI and the downstream
//! Sprinzl/modification analysis already consume.

use anyhow::{anyhow, Result};
use rayon::prelude::*;
use std::path::Path;

use ornament_alphabet::Dsq;
use ornament_scfg::model::nd;
use ornament_scfg::{
    align_glocal_banded, calc_qdb_bands, cm_pipeline_search, configure_local, configure_scores,
    parse_cm_file, parse_cm_records_file, Alignment, Cm, EmitMap, PipelineParams, QdbBands, Strand,
};

use super::stockholm::{AlignedRow, ModelMsa, ResCell};
use super::CMHit;

/// Bit-score floor handed to the scanner's overlap-resolution pass. Hits below this
/// are dropped before E-value filtering; it keeps significant tRNAs while discarding
/// the spurious low-scoring hits, mirroring the differential-test reporting floor.
const REPORTING_BITS: f32 = 20.0;

/// Per-model CM-DP memory budget (bytes). Large structured RNAs (e.g. LSU/SSU rRNA, `W` in the
/// thousands) blow up the scanning matrix — dominated by the BEGL rolling deck, which is
/// `O(W² · n_bifurcations)` — and would OOM a multi-model genome scan. Models whose estimate
/// exceeds this are skipped (with a warning) rather than crashing the whole run. ~1.5 GiB leaves
/// headroom under the default thread-count concurrency. (The proper fix for huge models is tighter
/// HMM bands / banded storage, not a bigger matrix.)
const MAX_SCAN_BYTES: usize = 1_500_000_000;

/// Estimate the peak CM scanning-DP footprint for `cm` at window `w_max` (the dominant `alpha` +
/// `alpha_begl` decks from `scan_core`). Used to skip pathologically large models before they OOM.
fn estimated_scan_bytes(cm: &ornament_scfg::Cm, w_max: usize) -> usize {
    let width = w_max + 1;
    let n_begl = cm
        .sttype
        .iter()
        .filter(|&&t| t == ornament_scfg::model::st::B)
        .count()
        .max(1);
    // alpha: 2 decks of m·width floats; alpha_begl: width·(n_begl·width) floats. f32 = 4 bytes.
    let cells = 2 * cm.m * width + width * width * n_begl;
    cells.saturating_mul(4)
}

/// Scan a FASTA file for CM hits using the native `ornament-scfg` engine.
///
/// `e_value` is the maximum E-value reported (mirrors `cmsearch -E`); hits from
/// calibrated models above it are filtered out. Uncalibrated models (no E-value)
/// fall back to the [`REPORTING_BITS`] floor only.
pub fn scan_native<P: AsRef<Path>, Q: AsRef<Path>>(
    cm_path: P,
    fasta: Q,
    e_value: f64,
) -> Result<Vec<CMHit>> {
    let cm_path = cm_path.as_ref();
    let fasta = fasta.as_ref();

    let mut cm = parse_cm_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM {}: {e}", cm_path.display()))?;
    // Glocal copy for per-hit alignment: `align_glocal` traces ROOT->END and needs an
    // un-localized model. The local model zeroes ROOT transitions (entry only via local
    // begins), so a glocal traceback over it would emit an all-insert parse (no match columns).
    let mut align_cm = cm.clone();
    configure_scores(&mut align_cm);
    configure_local(&mut cm); // cmsearch default (local mode)

    let records = ornament_alphabet::read_fasta(fasta)
        .map_err(|e| anyhow!("failed to read FASTA {}: {e}", fasta.display()))?;

    let w_max = cm.w as usize;
    let query_name = cm.name.clone();
    let query_acc = cm.acc.clone();

    // Emit map (node -> consensus column): built once, used to recover each hit's matched
    // consensus-column span (`mdl from`/`mdl to`) from its glocal traceback.
    let emap = EmitMap::build(&cm);

    // Query-dependent bands for the per-hit model-span alignment. The p7-filtered
    // `cm_pipeline_search` computes its own scan bands internally; `align_bands` (from the glocal
    // `align_cm`) bound the alignment over that glocal model — they must match the model the parse
    // runs over.
    let align_bands = calc_qdb_bands(&align_cm, QdbBands::DEFAULT_BETA);

    // Records are independent scans — search them in parallel. The per-record collect is
    // order-preserving, so the flattened result is identical regardless of thread count.
    let per_record: Vec<Vec<CMHit>> = records
        .par_iter()
        .map(|rec| -> Result<Vec<CMHit>> {
            let dsq = cm
                .abc
                .digitize(&rec.seq)
                .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;

            let (raw, _stats) =
                cm_pipeline_search(&cm, &dsq, w_max, REPORTING_BITS, PipelineParams::default())
                    .map_err(|e| anyhow!("pipeline failed for {}: {e}", cm.name))?;
            // Reverse-complement the record once if any minus-strand hit needs its model span.
            let rc = maybe_revcomp(&cm, &dsq, &raw, &rec.name)?;

            // Each hit's model-span recovery runs an independent banded glocal traceback; on a
            // long genome record with many hits this loop is the serial tail (the pipeline itself
            // is already parallel). Run the per-hit alignments in parallel — nested under the
            // per-record `par_iter`, rayon work-steals so a single big record still saturates
            // cores. `collect` preserves `raw`'s order, so output stays thread-count-independent.
            let rec_hits: Vec<CMHit> = raw
                .par_iter()
                .filter_map(|h| {
                    // Filter on E-value when the model is calibrated (mirrors `cmsearch -E`).
                    if let Some(ev) = h.evalue {
                        if ev > e_value {
                            return None;
                        }
                    }
                    let strand = match h.strand {
                        Strand::Plus => '+',
                        Strand::Minus => '-',
                    };
                    let (mdl_from, mdl_to) = hit_model_span(
                        &align_cm,
                        &dsq,
                        rc.as_deref(),
                        &emap,
                        &align_bands,
                        h.i,
                        h.j,
                        h.strand,
                    );
                    Some(CMHit {
                        target_name: rec.name.clone(),
                        target_start: h.i,
                        target_end: h.j,
                        strand,
                        query_name: query_name.clone(),
                        score: h.score as f64,
                        e_value: h.evalue.unwrap_or(f64::NAN),
                        gc_content: gc_fraction(&rec.seq, h.i, h.j),
                        mdl_from,
                        mdl_to,
                        query_accession: query_acc.clone(),
                        description: (!rec.desc.is_empty()).then(|| rec.desc.clone()),
                    })
                })
                .collect();
            Ok(rec_hits)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(per_record.into_iter().flatten().collect())
}

/// Scan a FASTA with **every** model in a (possibly multi-model) `.cm` file — the cmscan-style
/// "run the whole collection on a sequence" path. Work is parallelized across the full
/// `models × records` product, so a file of thousands of Rfam models keeps all cores busy even on
/// a single sequence (the across-model parallelism the per-model pipeline can't provide on its own).
///
/// `e_value` filters calibrated models' hits exactly as [`scan_native`]; each hit is tagged with
/// its model name (`CMHit::query_name`). Results are sorted best-first (E-value asc, then score).
pub fn scan_native_multi<P: AsRef<Path>, Q: AsRef<Path>>(
    cm_path: P,
    fasta: Q,
    e_value: f64,
) -> Result<Vec<CMHit>> {
    // Buffer every pair's hits, then sort best-first. Suited to interactive / small inputs.
    let collected = std::sync::Mutex::new(Vec::new());
    scan_multi_with(cm_path, fasta, e_value, |hits| {
        if !hits.is_empty() {
            collected.lock().unwrap().extend(hits);
        }
    })?;
    let mut hits = collected.into_inner().unwrap();
    // Best-first: E-value ascending (NaN/uncalibrated last), then score descending.
    hits.sort_by(|a: &CMHit, b: &CMHit| {
        let ea = if a.e_value.is_nan() {
            f64::INFINITY
        } else {
            a.e_value
        };
        let eb = if b.e_value.is_nan() {
            f64::INFINITY
        } else {
            b.e_value
        };
        ea.partial_cmp(&eb)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.score.total_cmp(&a.score))
    });
    Ok(hits)
}

/// Streaming variant of [`scan_native_multi`]: instead of buffering and sorting all hits, it hands
/// each (model, record) scan's hits to `on_hits` as that scan completes (arrival order, unsorted).
/// For long genome-scale runs over the whole Rfam this matters — results appear incrementally and
/// survive a crash, and the process never holds the full hit set. Returns the total hit count.
pub fn scan_native_multi_streaming<P, Q, F>(
    cm_path: P,
    fasta: Q,
    e_value: f64,
    on_hits: F,
) -> Result<usize>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    F: Fn(&[CMHit]) + Sync,
{
    use std::sync::atomic::{AtomicUsize, Ordering};
    let count = AtomicUsize::new(0);
    scan_multi_with(cm_path, fasta, e_value, |hits| {
        if !hits.is_empty() {
            count.fetch_add(hits.len(), Ordering::Relaxed);
            on_hits(&hits);
        }
    })?;
    Ok(count.load(Ordering::Relaxed))
}

/// Shared core: parse every model + memory-guard + per-model prep, then scan the feasible
/// model × record product in parallel, handing each pair's hits to `sink` as it completes. `sink`
/// is called once per scanned pair (possibly with an empty vec filtered by the callers) and must be
/// thread-safe.
fn scan_multi_with<P, Q, S>(cm_path: P, fasta: Q, e_value: f64, sink: S) -> Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    S: Fn(Vec<CMHit>) + Sync,
{
    let cm_path = cm_path.as_ref();
    let fasta = fasta.as_ref();

    let models = parse_cm_records_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM file {}: {e}", cm_path.display()))?;
    let records = ornament_alphabet::read_fasta(fasta)
        .map_err(|e| anyhow!("failed to read FASTA {}: {e}", fasta.display()))?;

    // Memory guard: skip models whose scanning DP would exceed the budget (huge rRNAs etc.) so one
    // pathological model can't OOM the whole collection scan. Reported, not silent.
    let mut n_skipped = 0usize;
    let feasible: Vec<bool> = models
        .iter()
        .map(|cm| {
            let ok = estimated_scan_bytes(cm, cm.w as usize) <= MAX_SCAN_BYTES;
            if !ok {
                n_skipped += 1;
                tracing::warn!(
                    "skipping model {} (W={}, ~{:.1} GiB DP > budget); too large for the scanning matrix",
                    cm.name,
                    cm.w,
                    estimated_scan_bytes(cm, cm.w as usize) as f64 / (1u64 << 30) as f64,
                );
            }
            ok
        })
        .collect();
    if n_skipped > 0 {
        tracing::warn!("{n_skipped} of {} models skipped (too large)", models.len());
    }

    // Per-model setup: a glocal-configured copy + emit map + QDB bands for the per-hit model-span
    // alignment. (The p7 pipeline computes its own scan bands internally; `align_bands` are the
    // glocal bands that bound the alignment over `align_cm`.)
    struct Prepared {
        cm: Cm,
        /// Glocal-configured copy used only for per-hit alignment (see `scan_native`).
        align_cm: Cm,
        align_bands: QdbBands,
        w_max: usize,
        emap: EmitMap,
    }
    let prepared: Vec<Prepared> = models
        .into_par_iter()
        .map(|mut cm| {
            let mut align_cm = cm.clone();
            configure_scores(&mut align_cm);
            configure_local(&mut cm); // cmsearch/cmscan default (local mode)
            let align_bands = calc_qdb_bands(&align_cm, QdbBands::DEFAULT_BETA);
            let w_max = cm.w as usize;
            let emap = EmitMap::build(&cm);
            Prepared {
                cm,
                align_cm,
                align_bands,
                w_max,
                emap,
            }
        })
        .collect();

    // Parallelize across the feasible model × record product: every (model, record) scan is
    // independent, so the whole collection saturates the thread pool regardless of how many
    // models or sequences there are. Each scan uses the p7-filtered pipeline (the brute QDB scan
    // is infeasible at genome scale), which prunes the genome cheaply before the CM DP.
    let pairs: Vec<(usize, usize)> = (0..prepared.len())
        .filter(|&mi| feasible[mi])
        .flat_map(|mi| (0..records.len()).map(move |ri| (mi, ri)))
        .collect();

    // Digitize (and reverse-complement) each record ONCE, up front. Both depend only on the
    // sequence and the alphabet — not the model — and every model parsed from a `.cm` file shares
    // the same RNA alphabet, so the digital forms are model-independent. This used to run inside
    // the per-(model × record) loop below: a genome scanned against thousands of Rfam models
    // re-digitized and re-reverse-complemented each sequence thousands of times. Owned `Vec`s (no
    // interior mutability) keep the downstream alignment DP reading plain slices. (Memory: holds
    // both strands of every record for the scan's duration — fine for the bacterial/chunked genomes
    // this path targets.)
    let digital: Vec<(Vec<Dsq>, Vec<Dsq>)> = if let Some(p0) = prepared.first() {
        let abc = &p0.cm.abc;
        records
            .par_iter()
            .map(|rec| -> Result<(Vec<Dsq>, Vec<Dsq>)> {
                let dsq = abc
                    .digitize(&rec.seq)
                    .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;
                let mut rc = dsq.clone();
                abc.revcomp(&mut rc)
                    .map_err(|e| anyhow!("failed to reverse-complement {}: {e}", rec.name))?;
                Ok((dsq, rc))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };

    pairs.par_iter().try_for_each(|&(mi, ri)| -> Result<()> {
        let p = &prepared[mi];
        let rec = &records[ri];
        let (dsq, rc) = &digital[ri];

        let (raw, _stats) = cm_pipeline_search(
            &p.cm,
            dsq,
            p.w_max,
            REPORTING_BITS,
            PipelineParams::default(),
        )
        .map_err(|e| anyhow!("pipeline failed for model {}: {e}", p.cm.name))?;

        let mut hits = Vec::new();
        for h in raw {
            if let Some(ev) = h.evalue {
                if ev > e_value {
                    continue;
                }
            }
            let strand = match h.strand {
                Strand::Plus => '+',
                Strand::Minus => '-',
            };
            let (mdl_from, mdl_to) = hit_model_span(
                &p.align_cm,
                dsq,
                Some(rc.as_slice()),
                &p.emap,
                &p.align_bands,
                h.i,
                h.j,
                h.strand,
            );
            hits.push(CMHit {
                target_name: rec.name.clone(),
                target_start: h.i,
                target_end: h.j,
                strand,
                query_name: p.cm.name.clone(),
                score: h.score as f64,
                e_value: h.evalue.unwrap_or(f64::NAN),
                gc_content: gc_fraction(&rec.seq, h.i, h.j),
                mdl_from,
                mdl_to,
                query_accession: p.cm.acc.clone(),
                description: (!rec.desc.is_empty()).then(|| rec.desc.clone()),
            });
        }
        sink(hits);
        Ok(())
    })
}

/// Scan a FASTA and return, per model, a Stockholm-ready multiple alignment of that model's
/// hits — the native equivalent of `cmsearch -A`. Like [`scan_native_multi`] it runs every
/// model in the `.cm` over every record (p7-filtered pipeline, same E-value/memory handling),
/// but keeps each surviving hit's full glocal parse instead of just its model span. Models with
/// no surviving hits are omitted, so the result is one [`ModelMsa`] (one Stockholm block) per
/// model that matched. Rows within a model are ordered best-first (E-value asc, then score).
pub fn scan_native_aligned<P: AsRef<Path>, Q: AsRef<Path>>(
    cm_path: P,
    fasta: Q,
    e_value: f64,
) -> Result<Vec<ModelMsa>> {
    let cm_path = cm_path.as_ref();
    let fasta = fasta.as_ref();

    let models = parse_cm_records_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM file {}: {e}", cm_path.display()))?;
    let records = ornament_alphabet::read_fasta(fasta)
        .map_err(|e| anyhow!("failed to read FASTA {}: {e}", fasta.display()))?;

    // Memory guard: skip models whose scanning DP would exceed the budget (see scan_native_multi).
    let feasible: Vec<bool> = models
        .iter()
        .map(|cm| {
            let ok = estimated_scan_bytes(cm, cm.w as usize) <= MAX_SCAN_BYTES;
            if !ok {
                tracing::warn!(
                    "skipping model {} (W={}, ~{:.1} GiB DP > budget); too large for the scanning matrix",
                    cm.name,
                    cm.w,
                    estimated_scan_bytes(cm, cm.w as usize) as f64 / (1u64 << 30) as f64,
                );
            }
            ok
        })
        .collect();

    struct Prepared {
        cm: Cm,
        align_cm: Cm,
        align_bands: QdbBands,
        w_max: usize,
        emap: EmitMap,
        rf: Vec<char>,
        ss_cons: Vec<char>,
    }
    let prepared: Vec<Prepared> = models
        .into_par_iter()
        .map(|mut cm| {
            let mut align_cm = cm.clone();
            configure_scores(&mut align_cm);
            configure_local(&mut cm);
            let align_bands = calc_qdb_bands(&align_cm, QdbBands::DEFAULT_BETA);
            let w_max = cm.w as usize;
            let emap = EmitMap::build(&cm);
            let (rf, ss_cons) = consensus_annotation(&cm, &emap);
            Prepared {
                cm,
                align_cm,
                align_bands,
                w_max,
                emap,
                rf,
                ss_cons,
            }
        })
        .collect();

    let pairs: Vec<(usize, usize)> = (0..prepared.len())
        .filter(|&mi| feasible[mi])
        .flat_map(|mi| (0..records.len()).map(move |ri| (mi, ri)))
        .collect();

    // Digitize + reverse-complement each record once, up front (see [`scan_multi_with`]): both are
    // model-independent, so doing them per (model × record) pair re-computed each record's digital
    // form once per model.
    let digital: Vec<(Vec<Dsq>, Vec<Dsq>)> = if let Some(p0) = prepared.first() {
        let abc = &p0.cm.abc;
        records
            .par_iter()
            .map(|rec| -> Result<(Vec<Dsq>, Vec<Dsq>)> {
                let dsq = abc
                    .digitize(&rec.seq)
                    .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;
                let mut rc = dsq.clone();
                abc.revcomp(&mut rc)
                    .map_err(|e| anyhow!("failed to reverse-complement {}: {e}", rec.name))?;
                Ok((dsq, rc))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };

    // Each entry: (model index, rows with their (E-value, -score) sort keys).
    type Keyed = (f64, f64, AlignedRow);
    let per_pair: Vec<(usize, Vec<Keyed>)> = pairs
        .par_iter()
        .map(|&(mi, ri)| -> Result<(usize, Vec<Keyed>)> {
            let p = &prepared[mi];
            let rec = &records[ri];
            let (dsq, rc) = &digital[ri];

            let (raw, _stats) = cm_pipeline_search(
                &p.cm,
                dsq,
                p.w_max,
                REPORTING_BITS,
                PipelineParams::default(),
            )
            .map_err(|e| anyhow!("pipeline failed for model {}: {e}", p.cm.name))?;

            let mut rows = Vec::new();
            for h in raw {
                if let Some(ev) = h.evalue {
                    if ev > e_value {
                        continue;
                    }
                }
                let (aln, seq) = align_hit_window(
                    &p.align_cm,
                    dsq,
                    Some(rc.as_slice()),
                    &p.emap,
                    &p.align_bands,
                    h.i,
                    h.j,
                    h.strand,
                );
                let label = format!("{}/{}-{}", rec.name, h.i, h.j);
                let row = aligned_row(&p.align_cm, &aln, seq, label);
                let evkey = h.evalue.unwrap_or(f64::INFINITY);
                rows.push((evkey, -(h.score as f64), row));
            }
            Ok((mi, rows))
        })
        .collect::<Result<Vec<_>>>()?;

    // Group rows back under their model.
    let mut by_model: Vec<Vec<Keyed>> = (0..prepared.len()).map(|_| Vec::new()).collect();
    for (mi, rows) in per_pair {
        by_model[mi].extend(rows);
    }

    let mut out = Vec::new();
    for (mi, mut rows) in by_model.into_iter().enumerate() {
        if rows.is_empty() {
            continue;
        }
        // Best-first: E-value ascending (uncalibrated last), then score descending.
        rows.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.1.total_cmp(&b.1))
        });
        let p = &prepared[mi];
        out.push(ModelMsa {
            model_name: p.cm.name.clone(),
            clen: p.cm.clen,
            rf: p.rf.clone(),
            ss_cons: p.ss_cons.clone(),
            rows: rows.into_iter().map(|(_, _, r)| r).collect(),
        });
    }
    Ok(out)
}

/// Align every sequence in `fasta` to the single model in `cm_path` over its full length — the
/// native `cmalign` analog. Unlike the scan paths this performs *no search*: each input record is
/// assumed to be a member of the model and aligned ROOT->END on the plus strand (glocal,
/// QDB-banded, like the per-hit alignment but spanning the whole sequence `1..=L`). Returns one
/// [`ModelMsa`] — a single Stockholm block with one row per input sequence, in input order, plus
/// the model's `#=GC RF` / `#=GC SS_cons` annotation.
///
/// Takes the first model only (`cmalign` is single-model); extra models in a collection `.cm`
/// are ignored.
pub fn align_sequences<P: AsRef<Path>, Q: AsRef<Path>>(cm_path: P, fasta: Q) -> Result<ModelMsa> {
    let cm_path = cm_path.as_ref();
    let fasta = fasta.as_ref();

    let mut cm = parse_cm_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM {}: {e}", cm_path.display()))?;
    // Glocal scores only (no local config): a full-length member aligns ROOT->END, the same
    // configuration `scan_native` uses for its per-hit alignment copy.
    configure_scores(&mut cm);
    let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);
    let emap = EmitMap::build(&cm);
    let (rf, ss_cons) = consensus_annotation(&cm, &emap);

    let records = ornament_alphabet::read_fasta(fasta)
        .map_err(|e| anyhow!("failed to read FASTA {}: {e}", fasta.display()))?;
    if records.is_empty() {
        return Err(anyhow!("no sequences to align in {}", fasta.display()));
    }

    // Each sequence aligns independently; saturate the pool over the input set.
    let rows = records
        .par_iter()
        .map(|rec| -> Result<AlignedRow> {
            let dsq = cm
                .abc
                .digitize(&rec.seq)
                .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;
            let l = dsq.len() - 2; // strip the index-0 / index-L+1 sentinels
            let aln = align_glocal_banded(&cm, &dsq, 1, l, &emap, &bands);
            let label = format!("{}/{}-{}", rec.name, 1, l);
            Ok(aligned_row(&cm, &aln, &dsq, label))
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(ModelMsa {
        model_name: cm.name.clone(),
        clen: cm.clen,
        rf,
        ss_cons,
        rows,
    })
}

/// Reverse-complement the (sentinel-padded) digital sequence once, but only if `raw` contains a
/// minus-strand hit that will need it for its model-span alignment. Returns `None` when no minus
/// hit is present (the common case) so plus-only records skip the allocation.
fn maybe_revcomp(
    cm: &Cm,
    dsq: &[Dsq],
    raw: &[ornament_scfg::Hit],
    name: &str,
) -> Result<Option<Vec<Dsq>>> {
    if !raw.iter().any(|h| matches!(h.strand, Strand::Minus)) {
        return Ok(None);
    }
    let mut v = dsq.to_vec();
    cm.abc
        .revcomp(&mut v)
        .map_err(|e| anyhow!("failed to reverse-complement {name}: {e}"))?;
    Ok(Some(v))
}

/// Matched consensus-column span (`mdl from`, `mdl to`) of a hit, recovered from its glocal
/// traceback. Plus hits align the original `dsq` over window `[i, j]`; minus hits (forward coords
/// `i > j`) align the reverse-complement `rc` over the mirrored window `[L-i+1, L-j+1]`, where the
/// Cap on the (QDB-banded) model-span alignment footprint. The banded matrix is
/// `(window+1)·Σ_v band(v)·12 B`; a full-length rRNA hit is still a few GB, and many such hits
/// align concurrently in a genome scan. Above the cap we report the coarse full-model span instead
/// — for a strong full-length rRNA-scale hit that is `1..clen` anyway, which the alignment confirms.
const MAX_SPAN_ALIGN_BYTES: u64 = 512 * 1024 * 1024;

/// Matched consensus-column span (`mdl from`, `mdl to`) of a hit, recovered from its glocal
/// traceback. Plus hits align the original `dsq` over window `[i, j]`; minus hits (forward coords
/// `i > j`) align the reverse-complement `rc` over the mirrored window `[L-i+1, L-j+1]`, where the
/// model columns already come out in 5'->3' model orientation. Returns `(0, 0)` if the parse
/// emitted no match columns. The alignment is QDB-banded (`O(m·window·band)`); for a hit whose
/// banded matrix would still exceed [`MAX_SPAN_ALIGN_BYTES`] it reports the coarse `(1, clen)` span.
#[allow(clippy::too_many_arguments)]
fn hit_model_span(
    cm: &Cm,
    dsq: &[Dsq],
    rc: Option<&[Dsq]>,
    emap: &EmitMap,
    bands: &QdbBands,
    i: usize,
    j: usize,
    strand: Strand,
) -> (usize, usize) {
    let window = i.abs_diff(j) as u64 + 1;
    let band_cells: u64 = (0..cm.m)
        .map(|v| (bands.dmax[v].saturating_sub(bands.dmin[v]) + 1) as u64)
        .sum();
    if (window + 1).saturating_mul(band_cells).saturating_mul(12) > MAX_SPAN_ALIGN_BYTES {
        return (1, cm.clen);
    }

    let (aln, _seq) = align_hit_window(cm, dsq, rc, emap, bands, i, j, strand);
    aln.residues
        .iter()
        .filter_map(|r| r.consensus)
        .fold(None, |acc: Option<(usize, usize)>, c| match acc {
            Some((lo, hi)) => Some((lo.min(c), hi.max(c))),
            None => Some((c, c)),
        })
        .unwrap_or((0, 0))
}

/// QDB-banded glocal alignment of a hit window, returning the parse and the digital sequence it was
/// aligned over. Plus hits align `dsq` over `[i, j]`; minus hits (forward coords `i > j`) align the
/// reverse-complement `rc` over the mirrored window, so the model columns come out in 5'->3'
/// model orientation either way. Shared by [`hit_model_span`] and the Stockholm output. The bands
/// bound memory to `O(m·window·band)` (vs the dense `O(m·window²)`) while giving an identical parse.
#[allow(clippy::too_many_arguments)]
fn align_hit_window<'a>(
    cm: &Cm,
    dsq: &'a [Dsq],
    rc: Option<&'a [Dsq]>,
    emap: &EmitMap,
    bands: &QdbBands,
    i: usize,
    j: usize,
    strand: Strand,
) -> (Alignment, &'a [Dsq]) {
    let l = dsq.len() - 2; // sequence length L (dsq has sentinels at index 0 and L+1)
    match strand {
        Strand::Plus => (align_glocal_banded(cm, dsq, i, j, emap, bands), dsq),
        Strand::Minus => {
            let rc = rc.expect("minus-strand hit requires a reverse-complemented sequence");
            (
                align_glocal_banded(cm, rc, l - i + 1, l - j + 1, emap, bands),
                rc,
            )
        }
    }
}

/// Build a Stockholm alignment row from a hit's glocal parse: each residue decoded to a base
/// character (uppercase for a match emission, lowercase for an insert) and tagged with its
/// consensus column or insert gap. `seq` is the digital sequence the parse was aligned over.
fn aligned_row(cm: &Cm, aln: &Alignment, seq: &[Dsq], label: String) -> AlignedRow {
    let residues = aln
        .residues
        .iter()
        .map(|r| {
            let raw = cm
                .abc
                .textize(&[seq[r.seq_pos]])
                .chars()
                .next()
                .unwrap_or('?');
            let (ch, consensus) = match r.consensus {
                Some(c) => (raw.to_ascii_uppercase(), Some(c)),
                None => (raw.to_ascii_lowercase(), None),
            };
            ResCell {
                consensus,
                insert_after: r.insert_after.unwrap_or(0),
                insert_right: r.insert_right,
                ch,
            }
        })
        .collect();
    AlignedRow { label, residues }
}

/// Per-consensus-column `#=GC RF` and `#=GC SS_cons` lines for a model (length `clen` each).
///
/// RF prefers the model's reference line, then its consensus string, then the consensus
/// residues assembled from the per-node annotation (`lcons`/`rcons` placed by the emit map —
/// case preserved, since `.cm` uses upper/lowercase to flag conservation), then a placeholder.
/// SS_cons is full WUSS (a port of Infernal's `CreateCMConsensus`): each MATP pair's bracket
/// family is set by its multifurcation order (`<>`/`()`/`[]`/`{}`) and unpaired MATL/MATR
/// columns are classified by their faces into `:` external, `_` hairpin, `-` bulge/interior,
/// or `,` multiloop.
fn consensus_annotation(cm: &Cm, emap: &EmitMap) -> (Vec<char>, Vec<char>) {
    let clen = cm.clen;
    let set = |v: &mut [char], col: i32, c: char| {
        let col = col as usize;
        if (1..=clen).contains(&col) {
            v[col - 1] = c;
        }
    };

    let rf: Vec<char> = if let Some(s) = cm.rf.as_ref().or(cm.consensus.as_ref()) {
        let mut v: Vec<char> = s.chars().collect();
        v.resize(clen, 'x');
        v
    } else {
        let mut v = vec!['x'; clen];
        for node in 0..cm.nodes {
            let ann = &cm.node_annot[node];
            match cm.ndtype[node] {
                nd::MATP => {
                    if let Some(c) = ann.lcons {
                        set(&mut v, emap.lpos[node], c);
                    }
                    if let Some(c) = ann.rcons {
                        set(&mut v, emap.rpos[node], c);
                    }
                }
                nd::MATL => {
                    if let Some(c) = ann.lcons {
                        set(&mut v, emap.lpos[node], c);
                    }
                }
                nd::MATR => {
                    if let Some(c) = ann.rcons {
                        set(&mut v, emap.rpos[node], c);
                    }
                }
                _ => {}
            }
        }
        v
    };

    let height = multifurcation_order(cm);
    let (inface, outface) = face_charts(cm);
    let mut ss = vec![':'; clen];
    for node in 0..cm.nodes {
        match cm.ndtype[node] {
            nd::MATP => {
                let (l, r) = match height[node] {
                    0 => ('<', '>'),
                    1 => ('(', ')'),
                    2 => ('[', ']'),
                    _ => ('{', '}'),
                };
                set(&mut ss, emap.lpos[node], l);
                set(&mut ss, emap.rpos[node], r);
            }
            nd::MATL => set(
                &mut ss,
                emap.lpos[node],
                face_char(inface[node], outface[node]),
            ),
            nd::MATR => set(
                &mut ss,
                emap.rpos[node],
                face_char(inface[node], outface[node]),
            ),
            _ => {}
        }
    }
    (rf, ss)
}

/// WUSS character for an unpaired (MATL/MATR) consensus column from its face counts
/// (Infernal `createFaceCharts`): external single-strand, hairpin loop, bulge/interior, or
/// multiloop.
fn face_char(inface: i32, outface: i32) -> char {
    if outface == 0 {
        ':' // external single-strand
    } else if inface == 0 && outface == 1 {
        '_' // hairpin loop
    } else if inface == 1 && outface == 1 {
        '-' // bulge / interior loop
    } else {
        ',' // multiloop
    }
}

/// Multifurcation order ("height") of the subtree at each node — a port of Infernal's
/// `createMultifurcationOrderChart`. Terminal stems are 0; a stem above a multifurcation into
/// terminal stems is 1; and so on. Selects the bracket family for each consensus pair.
fn multifurcation_order(cm: &Cm) -> Vec<i32> {
    let n = cm.nodes;
    let mut height = vec![0i32; n];
    let mut seg_has_pairs = vec![false; n];
    for node in (0..n).rev() {
        let v = cm.nodemap[node];
        seg_has_pairs[node] = match cm.ndtype[node] {
            nd::MATP => true,
            nd::END | nd::BIF => false,
            _ => seg_has_pairs[node + 1],
        };
        height[node] = match cm.ndtype[node] {
            nd::END => 0,
            nd::BIF => {
                let left = cm.ndidx[cm.cfirst[v] as usize];
                let right = cm.ndidx[cm.cnum[v] as usize];
                (height[left] + seg_has_pairs[left] as i32)
                    .max(height[right] + seg_has_pairs[right] as i32)
            }
            _ => height[node + 1],
        };
    }
    height
}

/// `inface`/`outface` face counts for every node — a port of Infernal's `createFaceCharts`.
/// `inface` counts faces in descendant subtrees (exclusive of the current pair); `outface`
/// counts faces above the node. Together they classify each unpaired column's loop context.
fn face_charts(cm: &Cm) -> (Vec<i32>, Vec<i32>) {
    let n = cm.nodes;
    let mut inface = vec![0i32; n];
    let mut outface = vec![0i32; n];

    for node in (0..n).rev() {
        let v = cm.nodemap[node];
        inface[node] = match cm.ndtype[node] {
            nd::END => 0,
            nd::BIF => {
                let left = cm.ndidx[cm.cfirst[v] as usize];
                let right = cm.ndidx[cm.cnum[v] as usize];
                inface[left] + inface[right]
            }
            _ if cm.ndtype[node + 1] == nd::MATP => 1,
            _ => inface[node + 1],
        };
    }

    for node in 0..n {
        let v = cm.nodemap[node];
        outface[node] = match cm.ndtype[node] {
            nd::ROOT => 0,
            nd::BEGL => {
                let parent = cm.ndidx[cm.plast[v] as usize];
                let y = cm.nodemap[parent];
                let right = cm.ndidx[cm.cnum[y] as usize];
                outface[parent] + inface[right]
            }
            nd::BEGR => {
                let parent = cm.ndidx[cm.plast[v] as usize];
                let y = cm.nodemap[parent];
                let left = cm.ndidx[cm.cfirst[y] as usize];
                outface[parent] + inface[left]
            }
            _ => {
                let parent = node - 1;
                if cm.ndtype[parent] == nd::MATP {
                    1
                } else {
                    outface[parent]
                }
            }
        };
    }

    (inface, outface)
}

/// GC fraction over the inclusive 1-based span `[min(i,j), max(i,j)]` of `seq`.
fn gc_fraction(seq: &str, i: usize, j: usize) -> f64 {
    let (lo, hi) = if i <= j { (i, j) } else { (j, i) };
    let bytes = seq.as_bytes();
    let start = lo.saturating_sub(1);
    let end = hi.min(bytes.len());
    if start >= end {
        return 0.0;
    }
    let span = &bytes[start..end];
    let gc = span
        .iter()
        .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / span.len() as f64
}
