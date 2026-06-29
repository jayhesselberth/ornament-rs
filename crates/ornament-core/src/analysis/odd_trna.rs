//! Odd tRNA detection

use super::{ModCompatibilityResult, TRNAHit};
use crate::modification::ModificationDatabase;

/// Detect odd tRNAs from a set of hits
pub fn detect_odd_trnas(
    hits: &[TRNAHit],
    mod_db: &ModificationDatabase,
    threshold: f64,
) -> Vec<ModCompatibilityResult> {
    hits.iter()
        .map(|hit| super::analyze_compatibility(hit, mod_db))
        .filter(|result| result.compatibility_score < threshold)
        .collect()
}
