//! modkit BedMethyl format parser
//!
//! Parses modification calls from modkit pileup output.

use serde::{Deserialize, Serialize};

/// A record from modkit BedMethyl output
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BedMethylRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub mod_code: String,
    pub score: u32,
    pub strand: char,
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
