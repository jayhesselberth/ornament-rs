//! Output format implementations

use crate::analysis::ModCompatibilityResult;

/// Convert results to JSON
pub fn to_json(results: &[ModCompatibilityResult]) -> serde_json::Result<String> {
    serde_json::to_string_pretty(results)
}

/// Convert results to TSV
pub fn to_tsv(results: &[ModCompatibilityResult]) -> String {
    let mut output = String::new();

    // Header
    output.push_str("id\tseq_name\tstart\tend\tstrand\tscore\tisotype\tanticodon\tis_odd\tcompatibility_score\n");

    for result in results {
        let hit = &result.hit;
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:?}\t{:.2}\t{}\t{}\t{}\t{:.3}\n",
            hit.id,
            hit.seq_name,
            hit.start,
            hit.end,
            hit.strand,
            hit.score,
            hit.isotype.as_deref().unwrap_or("-"),
            hit.anticodon.as_deref().unwrap_or("-"),
            result.is_odd,
            result.compatibility_score,
        ));
    }

    output
}
