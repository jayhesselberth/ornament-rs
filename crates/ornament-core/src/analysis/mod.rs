//! Analysis module
//!
//! Modification compatibility analysis and odd tRNA detection.

pub mod compatibility;
pub mod odd_trna;

use crate::SprinzlPosition;
use serde::{Deserialize, Serialize};

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
    /// Measured modifications (e.g. from modkit) joined onto Sprinzl positions, with a
    /// verdict per call. Empty when no modkit data was supplied.
    #[serde(default)]
    pub measured: Vec<MeasuredModification>,
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

/// A measured modification (e.g. a modkit call) joined onto a Sprinzl position, together
/// with how it squares against the genomic sequence and the MODOMICS expectations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeasuredModification {
    pub position: SprinzlPosition,
    /// modkit modification code as reported in the bedMethyl.
    pub mod_code: String,
    /// Genomic base underlying the call at this position.
    pub observed_base: crate::RnaBase,
    pub mod_frequency: f64,
    pub coverage: u32,
    pub verdict: MeasuredVerdict,
}

/// Verdict for a measured modification relative to sequence + MODOMICS expectations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeasuredVerdict {
    /// Chemically possible on this base and expected by MODOMICS at this position.
    Consistent,
    /// Chemically possible, but MODOMICS does not expect this modification here.
    Unexpected,
    /// The underlying base cannot carry this modification (impossible chemistry).
    Incompatible,
}

pub use compatibility::{
    analyze_batch, analyze_compatibility, analyze_compatibility_with_mods, BatchAnalysisResult,
};
pub use odd_trna::detect_odd_trnas;
