//! Infernal output parser
//!
//! Parses cmsearch tabular and Stockholm alignment outputs.

use serde::{Deserialize, Serialize};

/// A covariance model hit from cmsearch
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CMHit {
    pub target_name: String,
    pub target_start: usize,
    pub target_end: usize,
    pub strand: char,
    pub query_name: String,
    pub score: f64,
    pub e_value: f64,
    pub gc_content: f64,
}

/// Alignment from cmsearch Stockholm output
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CMAlignment {
    pub hit: CMHit,
    pub target_seq: String,
    pub consensus_seq: String,
    pub structure: String,
}

/// Parse cmsearch tabular output (--tblout)
pub fn parse_tblout(content: &str) -> Vec<CMHit> {
    let mut hits = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 18 {
            continue;
        }

        // tblout format: target name, accession, query name, accession, mdl, mdl from, mdl to,
        // seq from, seq to, strand, trunc, pass, gc, bias, score, E-value, inc, description
        if let (Ok(start), Ok(end), Ok(score), Ok(e_value), Ok(gc)) = (
            fields[7].parse::<usize>(),
            fields[8].parse::<usize>(),
            fields[14].parse::<f64>(),
            fields[15].parse::<f64>(),
            fields[12].parse::<f64>(),
        ) {
            hits.push(CMHit {
                target_name: fields[0].to_string(),
                target_start: start,
                target_end: end,
                strand: fields[9].chars().next().unwrap_or('+'),
                query_name: fields[2].to_string(),
                score,
                e_value,
                gc_content: gc,
            });
        }
    }

    hits
}
