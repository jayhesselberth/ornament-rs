//! Modification compatibility analysis

use super::{
    MeasuredModification, MeasuredVerdict, ModCompatibilityResult, ModificationIncompatibility,
    Severity, TRNAHit,
};
use crate::integration::modkit::ObservedModCall;
use crate::modification::Isotype;
use crate::modification::{ModificationDatabase, SprinzlMapper};
use crate::{ConservationLevel, RnaBase, SprinzlPosition};
use std::collections::HashMap;

/// Analyze modification compatibility for a tRNA hit (sequence only, no measured data).
pub fn analyze_compatibility(
    hit: &TRNAHit,
    mod_db: &ModificationDatabase,
) -> ModCompatibilityResult {
    analyze_compatibility_with_mods(hit, mod_db, &HashMap::new())
}

/// Analyze modification compatibility, additionally folding in measured modifications
/// (e.g. modkit calls joined onto Sprinzl positions via `integration::join_to_sprinzl`).
///
/// Behaves identically to [`analyze_compatibility`] when `observed` is empty. Each measured
/// call is graded against the genomic base (impossible chemistry -> `Incompatible`) and the
/// MODOMICS expectation (`Consistent` vs `Unexpected`). A call that is `Incompatible` or
/// `Unexpected` marks the tRNA "odd", in addition to the sequence-based criteria.
pub fn analyze_compatibility_with_mods(
    hit: &TRNAHit,
    mod_db: &ModificationDatabase,
    observed: &HashMap<SprinzlPosition, ObservedModCall>,
) -> ModCompatibilityResult {
    let mapper = SprinzlMapper::new_standard();

    // Map the sequence to Sprinzl positions using the structure as alignment guide
    // The structure string from cmsearch corresponds to CM columns
    let sprinzl_alignment = map_sequence_to_sprinzl(hit, &mapper);

    // Check each position for modification compatibility
    let mut incompatibilities = Vec::new();
    let mut positions_checked = 0;
    let mut positions_compatible = 0;

    // Get isotype for isotype-specific checks
    let isotype = hit.isotype.as_ref().map(Isotype::new);

    for (sprinzl_pos, seq_idx) in &sprinzl_alignment {
        // Get the base at this position
        let base_char = hit.sequence.chars().nth(*seq_idx);
        let observed_base = base_char.and_then(RnaBase::from_dna_char);

        let Some(observed) = observed_base else {
            continue; // Skip non-standard bases (N, etc.)
        };

        // Get expected modifications at this position
        let expectations = if let Some(ref iso) = isotype {
            mod_db.get_expectations_for_isotype(sprinzl_pos, iso)
        } else {
            mod_db.get_expectations(sprinzl_pos)
        };

        if expectations.is_empty() {
            continue; // No modification expected here
        }

        positions_checked += 1;

        // Check if the observed base is compatible with any expected modification
        let mut position_compatible = false;

        for expectation in &expectations {
            for modification in &expectation.modifications {
                if modification.is_compatible(observed) {
                    position_compatible = true;
                    break;
                }

                // Found an incompatibility
                let severity = match expectation.conservation {
                    ConservationLevel::Universal => Severity::Critical,
                    ConservationLevel::DomainSpecific => Severity::Major,
                    ConservationLevel::IsotypeSpecific => Severity::Major,
                    ConservationLevel::Rare => Severity::Minor,
                };

                incompatibilities.push(ModificationIncompatibility {
                    position: sprinzl_pos.clone(),
                    observed_base: observed,
                    expected_mod_name: modification.short_name.clone(),
                    severity,
                });
            }

            if position_compatible {
                break;
            }
        }

        if position_compatible {
            positions_compatible += 1;
        }
    }

    // Calculate compatibility score
    let compatibility_score = if positions_checked > 0 {
        positions_compatible as f64 / positions_checked as f64
    } else {
        1.0 // No positions to check = fully compatible
    };

    // Determine if this is an "odd" tRNA
    // Odd if: score < 1.0 AND has critical/major incompatibilities
    let has_significant_incompatibility = incompatibilities
        .iter()
        .any(|i| matches!(i.severity, Severity::Critical | Severity::Major));
    let mut is_odd = compatibility_score < 1.0 && has_significant_incompatibility;

    // Fold in measured modifications (e.g. modkit calls). A call that is impossible on the
    // underlying base, or unexpected vs MODOMICS, also makes the tRNA "odd".
    let measured = grade_measured_mods(hit, mod_db, &sprinzl_alignment, isotype.as_ref(), observed);
    if measured
        .iter()
        .any(|m| !matches!(m.verdict, MeasuredVerdict::Consistent))
    {
        is_odd = true;
    }

    ModCompatibilityResult {
        hit: hit.clone(),
        sprinzl_alignment,
        incompatibilities,
        is_odd,
        compatibility_score,
        measured,
    }
}

/// Grade each measured modification call against the genomic base and MODOMICS expectations.
fn grade_measured_mods(
    hit: &TRNAHit,
    mod_db: &ModificationDatabase,
    sprinzl_alignment: &HashMap<SprinzlPosition, usize>,
    isotype: Option<&Isotype>,
    observed: &HashMap<SprinzlPosition, ObservedModCall>,
) -> Vec<MeasuredModification> {
    let mut measured = Vec::new();
    for (pos, call) in observed {
        // Locate the genomic base under this measured call.
        let Some(&seq_idx) = sprinzl_alignment.get(pos) else {
            continue; // position not aligned in this hit
        };
        let Some(base) = hit
            .sequence
            .chars()
            .nth(seq_idx)
            .and_then(RnaBase::from_dna_char)
        else {
            continue;
        };

        let verdict = match mod_db.get_by_modkit_code(&call.mod_code) {
            // Impossible chemistry: the called modification can't sit on this base.
            Some(m) if !m.is_compatible(base) => MeasuredVerdict::Incompatible,
            // Possible: is it expected here by MODOMICS?
            Some(m) => {
                let expectations = if let Some(iso) = isotype {
                    mod_db.get_expectations_for_isotype(pos, iso)
                } else {
                    mod_db.get_expectations(pos)
                };
                let expected = expectations
                    .iter()
                    .any(|e| e.modifications.iter().any(|em| em.short_name == m.short_name));
                if expected {
                    MeasuredVerdict::Consistent
                } else {
                    MeasuredVerdict::Unexpected
                }
            }
            // Unknown modkit code: can't assess chemistry, so treat as unexpected.
            None => MeasuredVerdict::Unexpected,
        };

        measured.push(MeasuredModification {
            position: pos.clone(),
            mod_code: call.mod_code.clone(),
            observed_base: base,
            mod_frequency: call.mod_frequency,
            coverage: call.coverage,
            verdict,
        });
    }
    // Stable, position-sorted order keeps results deterministic across HashMap iteration.
    measured.sort_by(|a, b| a.position.0.cmp(&b.position.0));
    measured
}

/// Map a tRNA sequence to Sprinzl positions.
///
/// Primary path: align the hit sequence to the native Sprinzl reference covariance model
/// (`modification::align_to_sprinzl`), which gives a real consensus-column → Sprinzl
/// mapping that handles insertions/deletions in the D-loop and variable arm. Only if that
/// produces nothing (e.g. a non-tRNA / non-standard-residue sequence the model can't
/// digitize) do we fall back to the naive 1:1 mapping.
fn map_sequence_to_sprinzl(
    hit: &TRNAHit,
    mapper: &SprinzlMapper,
) -> HashMap<SprinzlPosition, usize> {
    let aligned = crate::modification::align_to_sprinzl(&hit.sequence);
    if !aligned.is_empty() {
        return aligned;
    }

    // Fallback: assume sequence is already aligned to standard positions (1:1, ungapped).
    let mut result = HashMap::new();
    for (seq_idx, _) in hit.sequence.chars().enumerate() {
        if let Some(sprinzl) = mapper.get_sprinzl(seq_idx) {
            result.insert(sprinzl.clone(), seq_idx);
        }
    }
    result
}

/// Analyze multiple tRNA hits and return summary statistics
pub fn analyze_batch(hits: &[TRNAHit], mod_db: &ModificationDatabase) -> BatchAnalysisResult {
    let results: Vec<ModCompatibilityResult> = hits
        .iter()
        .map(|hit| analyze_compatibility(hit, mod_db))
        .collect();

    let total = results.len();
    let odd_count = results.iter().filter(|r| r.is_odd).count();
    let avg_score = if total > 0 {
        results.iter().map(|r| r.compatibility_score).sum::<f64>() / total as f64
    } else {
        1.0
    };

    BatchAnalysisResult {
        results,
        total_trnas: total,
        odd_trnas: odd_count,
        average_compatibility: avg_score,
    }
}

/// Result of batch analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BatchAnalysisResult {
    pub results: Vec<ModCompatibilityResult>,
    pub total_trnas: usize,
    pub odd_trnas: usize,
    pub average_compatibility: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::Strand;

    #[test]
    fn test_analyze_compatibility_compatible() {
        let hit = TRNAHit {
            id: "test1".to_string(),
            seq_name: "chr1".to_string(),
            start: 1000,
            end: 1072,
            strand: Strand::Plus,
            score: 80.0,
            isotype: Some("Ala".to_string()),
            anticodon: Some("AGC".to_string()),
            // U at position 55 (compatible with Psi)
            sequence:
                "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
                    .to_string(),
            structure:
                "(((((((..((((.........)))).(((((.......))))).....(((((.......))))))))))))...."
                    .to_string(),
        };

        let db = ModificationDatabase::eukaryotic();
        let result = analyze_compatibility(&hit, &db);

        // Should have some mapping
        assert!(!result.sprinzl_alignment.is_empty());
    }

    #[test]
    fn test_analyze_batch() {
        let hits = vec![TRNAHit {
            id: "test1".to_string(),
            seq_name: "chr1".to_string(),
            start: 1000,
            end: 1072,
            strand: Strand::Plus,
            score: 80.0,
            isotype: Some("Ala".to_string()),
            anticodon: Some("AGC".to_string()),
            sequence:
                "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
                    .to_string(),
            structure: "".to_string(),
        }];

        let db = ModificationDatabase::eukaryotic();
        let batch_result = analyze_batch(&hits, &db);

        assert_eq!(batch_result.total_trnas, 1);
    }
}
