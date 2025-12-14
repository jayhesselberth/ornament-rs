//! MODOMICS database parser
//!
//! Parses modification data from the MODOMICS REST API JSON format.
//! Download data from: https://genesilico.pl/modomics/api/modifications

use super::types::{Modification, ModCode, RnaBase};
use rustc_hash::FxHashMap;
use serde::Deserialize;
use std::path::Path;

/// Raw MODOMICS modification entry from JSON
#[derive(Debug, Deserialize)]
pub struct ModomicsEntry {
    pub id: u32,
    pub name: String,
    pub short_name: String,
    #[serde(default)]
    pub new_abbrev: Option<String>,
    #[serde(default)]
    pub reference_moiety: Vec<String>,
    #[serde(default)]
    pub formula: Option<String>,
    #[serde(default)]
    pub mass_avg: Option<f64>,
    #[serde(default)]
    pub smile: Option<String>,
}

/// Parse MODOMICS JSON file into modification entries
pub fn parse_modomics_file(path: &Path) -> Result<FxHashMap<String, Modification>, ModomicsError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| ModomicsError::IoError(e.to_string()))?;
    parse_modomics_json(&content)
}

/// Parse MODOMICS JSON string into modification entries
pub fn parse_modomics_json(json: &str) -> Result<FxHashMap<String, Modification>, ModomicsError> {
    // MODOMICS returns a map with string keys (modification IDs)
    let raw: FxHashMap<String, ModomicsEntry> = serde_json::from_str(json)
        .map_err(|e| ModomicsError::ParseError(e.to_string()))?;

    let mut modifications = FxHashMap::default();

    for (_id, entry) in raw {
        if let Some(modification) = convert_entry(&entry) {
            modifications.insert(modification.short_name.clone(), modification);
        }
    }

    Ok(modifications)
}

/// Convert a MODOMICS entry to our Modification type
fn convert_entry(entry: &ModomicsEntry) -> Option<Modification> {
    // Determine parent base from reference_moiety
    let parent_base = match entry.reference_moiety.first().map(|s| s.as_str()) {
        Some("A") => RnaBase::A,
        Some("C") => RnaBase::C,
        Some("G") => RnaBase::G,
        Some("U") => RnaBase::U,
        _ => return None, // Skip entries without clear parent base (X, QtRNA, etc.)
    };

    // Determine incompatible bases (all bases except parent)
    let incompatible_bases: Vec<RnaBase> = [RnaBase::A, RnaBase::C, RnaBase::G, RnaBase::U]
        .into_iter()
        .filter(|&b| b != parent_base)
        .collect();

    // Get unicode character if available
    let modomics_unicode = entry.new_abbrev.as_ref()
        .and_then(|s| s.chars().next());

    // Create primary code
    let code = if let Some(unicode) = modomics_unicode {
        ModCode::Unicode(unicode)
    } else if entry.short_name.len() == 1 {
        ModCode::SingleChar(entry.short_name.chars().next().unwrap())
    } else {
        ModCode::ShortName(entry.short_name.clone())
    };

    // Build alternative codes
    let mut alt_codes = Vec::new();
    if entry.short_name.len() == 1 {
        if let Some(c) = entry.short_name.chars().next() {
            if !matches!(code, ModCode::SingleChar(_)) {
                alt_codes.push(ModCode::SingleChar(c));
            }
        }
    }

    Some(Modification {
        name: entry.name.clone(),
        short_name: entry.short_name.clone(),
        code,
        alt_codes,
        parent_base,
        genomic_expectation: parent_base,
        incompatible_bases,
        chebi_id: None, // MODOMICS doesn't include ChEBI in basic API
        modomics_unicode,
    })
}

/// Errors from MODOMICS parsing
#[derive(Debug)]
pub enum ModomicsError {
    IoError(String),
    ParseError(String),
}

impl std::fmt::Display for ModomicsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ModomicsError::IoError(e) => write!(f, "IO error: {}", e),
            ModomicsError::ParseError(e) => write!(f, "Parse error: {}", e),
        }
    }
}

impl std::error::Error for ModomicsError {}

/// Mapping from common names/aliases to MODOMICS short names
pub fn get_common_aliases() -> FxHashMap<&'static str, &'static str> {
    let mut aliases = FxHashMap::default();

    // Pseudouridine aliases
    aliases.insert("Psi", "Y");
    aliases.insert("psi", "Y");
    aliases.insert("Ψ", "Y");

    // Other common aliases
    aliases.insert("rT", "m5U");  // ribothymidine
    aliases.insert("T", "m5U");

    aliases
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_modomics_entry() {
        let json = r#"{
            "118": {
                "id": 118,
                "name": "pseudouridine",
                "short_name": "Y",
                "new_abbrev": "P",
                "reference_moiety": ["U"],
                "formula": "C9H12N2O6"
            }
        }"#;

        let result = parse_modomics_json(json).unwrap();
        assert!(result.contains_key("Y"));

        let psi = result.get("Y").unwrap();
        assert_eq!(psi.name, "pseudouridine");
        assert_eq!(psi.parent_base, RnaBase::U);
        assert!(psi.incompatible_bases.contains(&RnaBase::A));
        assert!(!psi.incompatible_bases.contains(&RnaBase::U));
    }

    #[test]
    fn test_convert_entry_filters_invalid() {
        // Entry with unknown reference base should be filtered
        let entry = ModomicsEntry {
            id: 999,
            name: "unknown".to_string(),
            short_name: "X".to_string(),
            new_abbrev: None,
            reference_moiety: vec!["X".to_string()],
            formula: None,
            mass_avg: None,
            smile: None,
        };

        assert!(convert_entry(&entry).is_none());
    }
}
