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
use ornament_scfg::{
    align_glocal, calc_qdb_bands, configure_local, configure_scores, cyk_search_banded,
    parse_cm_file, parse_cm_records_file, Cm, EmitMap, QdbBands, Strand,
};

use super::CMHit;

/// Bit-score floor handed to the scanner's overlap-resolution pass. Hits below this
/// are dropped before E-value filtering; it keeps significant tRNAs while discarding
/// the spurious low-scoring hits, mirroring the differential-test reporting floor.
const REPORTING_BITS: f32 = 20.0;

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

    // Query-dependent bands: computed once for the model and shared read-only across every
    // record scan. They give identical hits to the unbanded scan at a fraction of the cost.
    let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);

    // Records are independent scans — search them in parallel. The per-record collect is
    // order-preserving, so the flattened result is identical regardless of thread count.
    let per_record: Vec<Vec<CMHit>> = records
        .par_iter()
        .map(|rec| -> Result<Vec<CMHit>> {
            let dsq = cm
                .abc
                .digitize(&rec.seq)
                .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;

            let raw = cyk_search_banded(&cm, &dsq, w_max, REPORTING_BITS, &bands);
            // Reverse-complement the record once if any minus-strand hit needs its model span.
            let rc = maybe_revcomp(&cm, &dsq, &raw, &rec.name)?;

            let mut rec_hits = Vec::new();
            for h in raw {
                // Filter on E-value when the model is calibrated (mirrors `cmsearch -E`).
                if let Some(ev) = h.evalue {
                    if ev > e_value {
                        continue;
                    }
                }
                let strand = match h.strand {
                    Strand::Plus => '+',
                    Strand::Minus => '-',
                };
                let (mdl_from, mdl_to) =
                    hit_model_span(&align_cm, &dsq, rc.as_deref(), &emap, h.i, h.j, h.strand);
                rec_hits.push(CMHit {
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
                });
            }
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
    let cm_path = cm_path.as_ref();
    let fasta = fasta.as_ref();

    let models = parse_cm_records_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM file {}: {e}", cm_path.display()))?;
    let records = ornament_alphabet::read_fasta(fasta)
        .map_err(|e| anyhow!("failed to read FASTA {}: {e}", fasta.display()))?;

    // Per-model setup (QDB bands + max window + emit map), computed once per model in parallel.
    struct Prepared {
        cm: Cm,
        /// Glocal-configured copy used only for per-hit alignment (see `scan_native`).
        align_cm: Cm,
        bands: QdbBands,
        w_max: usize,
        emap: EmitMap,
    }
    let prepared: Vec<Prepared> = models
        .into_par_iter()
        .map(|mut cm| {
            let mut align_cm = cm.clone();
            configure_scores(&mut align_cm);
            configure_local(&mut cm); // cmsearch/cmscan default (local mode)
            let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);
            let w_max = cm.w as usize;
            let emap = EmitMap::build(&cm);
            Prepared {
                cm,
                align_cm,
                bands,
                w_max,
                emap,
            }
        })
        .collect();

    // Parallelize across the full model × record product: every (model, record) scan is
    // independent, so the whole collection saturates the thread pool regardless of how many
    // models or sequences there are. Order-preserving collect → thread-count-independent output.
    let pairs: Vec<(usize, usize)> = (0..prepared.len())
        .flat_map(|mi| (0..records.len()).map(move |ri| (mi, ri)))
        .collect();

    let per_pair: Vec<Vec<CMHit>> = pairs
        .par_iter()
        .map(|&(mi, ri)| -> Result<Vec<CMHit>> {
            let p = &prepared[mi];
            let rec = &records[ri];
            let dsq =
                p.cm.abc
                    .digitize(&rec.seq)
                    .map_err(|e| anyhow!("failed to digitize sequence {}: {e}", rec.name))?;

            let raw = cyk_search_banded(&p.cm, &dsq, p.w_max, REPORTING_BITS, &p.bands);
            let rc = maybe_revcomp(&p.cm, &dsq, &raw, &rec.name)?;

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
                    &dsq,
                    rc.as_deref(),
                    &p.emap,
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
            Ok(hits)
        })
        .collect::<Result<Vec<_>>>()?;

    let mut hits: Vec<CMHit> = per_pair.into_iter().flatten().collect();
    // Best-first: E-value ascending (NaN/uncalibrated last), then score descending.
    hits.sort_by(|a, b| {
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
/// model columns already come out in 5'->3' model orientation. Returns `(0, 0)` if the parse
/// emitted no match columns (defensive — real hits always match at least one).
fn hit_model_span(
    cm: &Cm,
    dsq: &[Dsq],
    rc: Option<&[Dsq]>,
    emap: &EmitMap,
    i: usize,
    j: usize,
    strand: Strand,
) -> (usize, usize) {
    let l = dsq.len() - 2; // sequence length L (dsq has sentinels at index 0 and L+1)
    let aln = match strand {
        Strand::Plus => align_glocal(cm, dsq, i, j, emap),
        Strand::Minus => {
            let rc = rc.expect("minus-strand hit requires a reverse-complemented sequence");
            align_glocal(cm, rc, l - i + 1, l - j + 1, emap)
        }
    };
    aln.residues
        .iter()
        .filter_map(|r| r.consensus)
        .fold(None, |acc: Option<(usize, usize)>, c| match acc {
            Some((lo, hi)) => Some((lo.min(c), hi.max(c))),
            None => Some((c, c)),
        })
        .unwrap_or((0, 0))
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
