//! Analysis module
//!
//! Modification compatibility analysis and odd tRNA detection.

pub mod compatibility;
pub mod odd_trna;

use serde::{Deserialize, Serialize};
use crate::SprinzlPosition;

/// Represents a tRNA hit with associated metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TRNAHit {
    pub id: String,
    pub seq_name: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub score: f64,
    pub isotype: Option<String>,
    pub anticodon: Option<String>,
    pub sequence: String,
    pub structure: String,
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl From<char> for Strand {
    fn from(c: char) -> Self {
        match c {
            '-' => Strand::Minus,
            _ => Strand::Plus,
        }
    }
}

/// Result of modification compatibility analysis for a tRNA
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModCompatibilityResult {
    pub hit: TRNAHit,
    pub sprinzl_alignment: std::collections::HashMap<SprinzlPosition, usize>,
    pub incompatibilities: Vec<ModificationIncompatibility>,
    pub is_odd: bool,
    pub compatibility_score: f64,
}

/// A specific modification incompatibility found at a position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModificationIncompatibility {
    pub position: SprinzlPosition,
    pub observed_base: crate::RnaBase,
    pub expected_mod_name: String,
    pub severity: Severity,
}

/// Severity of a modification incompatibility
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Severity {
    Critical,
    Major,
    Minor,
}

pub use compatibility::analyze_compatibility;
pub use odd_trna::detect_odd_trnas;
