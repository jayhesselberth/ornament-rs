//! modkit BedMethyl format parser
//!
//! Parses modification calls from modkit pileup output and joins them onto the
//! native Sprinzl positions of a tRNA hit (via `modification::align_to_sprinzl`).

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::analysis::{Strand, TRNAHit};
use crate::modification::{align_to_sprinzl, SprinzlPosition};

/// A record from modkit BedMethyl output.
///
/// modkit's (extended) bedMethyl has 18 columns; we keep the ones we need:
/// chrom, start, end, mod-code, score, strand, then (cols 10/11) Nvalid_cov and
/// the percent-modified fraction. `start` is 0-based per the BED convention.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BedMethylRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub mod_code: String,
    pub score: u32,
    pub strand: char,
    /// Nvalid_cov: reads with a valid call at this position.
    pub coverage: u32,
    /// Percent modified (0-100), as emitted by modkit.
    pub mod_frequency: f64,
}

/// A modkit call joined onto a single Sprinzl position of a tRNA.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservedModCall {
    /// modkit modification code as written in the bedMethyl (e.g. "m", "a", "17802").
    pub mod_code: String,
    pub coverage: u32,
    pub mod_frequency: f64,
}

/// Parse a BedMethyl file
pub fn parse_bedmethyl(content: &str) -> Vec<BedMethylRecord> {
    let mut records = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        if let (Ok(start), Ok(end), Ok(score), Ok(coverage), Ok(freq)) = (
            fields[1].parse::<usize>(),
            fields[2].parse::<usize>(),
            fields[4].parse::<u32>(),
            fields[9].parse::<u32>(),
            fields[10].parse::<f64>(),
        ) {
            records.push(BedMethylRecord {
                chrom: fields[0].to_string(),
                start,
                end,
                mod_code: fields[3].to_string(),
                score,
                strand: fields[5].chars().next().unwrap_or('+'),
                coverage,
                mod_frequency: freq,
            });
        }
    }

    records
}

/// Map a BedMethyl record's genomic position to a 0-based index into a hit's sequence.
///
/// `TRNAHit::start`/`end` follow the cmsearch convention: 1-based inclusive coordinates
/// with `start` marking the 5' residue (the larger coordinate on the minus strand).
/// BedMethyl `start` is 0-based, so its 1-based genomic position is `record.start + 1`.
fn record_seq_index(hit: &TRNAHit, record: &BedMethylRecord) -> Option<usize> {
    let genomic = record.start + 1; // 0-based BED -> 1-based genomic
    match hit.strand {
        Strand::Plus => genomic.checked_sub(hit.start),
        Strand::Minus => hit.start.checked_sub(genomic),
    }
}

/// Join modkit BedMethyl records onto the Sprinzl positions of a single tRNA hit.
///
/// Each record is matched to the hit by chromosome and strand, its genomic position is
/// translated to a position within the hit's sequence, and that position is mapped to a
/// Sprinzl label via the native reference alignment. Records that fall outside the hit or
/// land on insert residues (no consensus column) are dropped. If two records hit the same
/// Sprinzl position the higher-frequency call wins (keeps the join deterministic).
pub fn join_to_sprinzl(
    hit: &TRNAHit,
    records: &[BedMethylRecord],
) -> HashMap<SprinzlPosition, ObservedModCall> {
    // Sprinzl label -> seq index, inverted to seq index -> Sprinzl label.
    let sprinzl = align_to_sprinzl(&hit.sequence);
    let mut idx_to_pos: HashMap<usize, SprinzlPosition> = HashMap::new();
    for (pos, idx) in sprinzl {
        idx_to_pos.insert(idx, pos);
    }

    let mut out: HashMap<SprinzlPosition, ObservedModCall> = HashMap::new();
    for record in records {
        if record.chrom != hit.seq_name {
            continue;
        }
        if Strand::from(record.strand) != hit.strand {
            continue;
        }
        let Some(seq_idx) = record_seq_index(hit, record) else {
            continue;
        };
        let Some(pos) = idx_to_pos.get(&seq_idx) else {
            continue; // outside the aligned region or an insert residue
        };

        let call = ObservedModCall {
            mod_code: record.mod_code.clone(),
            coverage: record.coverage,
            mod_frequency: record.mod_frequency,
        };
        out.entry(pos.clone())
            .and_modify(|existing| {
                if call.mod_frequency > existing.mod_frequency {
                    *existing = call.clone();
                }
            })
            .or_insert(call);
    }
    out
}
