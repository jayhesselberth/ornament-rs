//! Sprinzl position mapping
//!
//! Maps covariance model alignment positions to standard Sprinzl tRNA positions.
//! Sprinzl numbering uses positions 1-76 with variable regions using letters (e.g., 17a).

use std::collections::HashMap;
use super::types::SprinzlPosition;

/// Maps CM alignment positions to Sprinzl positions
pub struct SprinzlMapper {
    /// Map from CM column index to Sprinzl position
    cm_to_sprinzl: HashMap<usize, SprinzlPosition>,
    /// Map from Sprinzl position back to CM column index
    sprinzl_to_cm: HashMap<SprinzlPosition, usize>,
}

impl SprinzlMapper {
    /// Create a new mapper with the standard tRNA Sprinzl mapping
    pub fn new_standard() -> Self {
        let mut cm_to_sprinzl = HashMap::new();
        let mut sprinzl_to_cm = HashMap::new();

        // Standard tRNA Sprinzl positions (71 core + variable loop)
        // Based on the canonical tRNA structure
        let positions = vec![
            // Acceptor stem (1-7, 66-72)
            "1", "2", "3", "4", "5", "6", "7",
            // D-stem/loop (8-16, 21-25)
            "8", "9", "10", "11", "12", "13", "14", "15", "16",
            // Variable D-loop positions
            "17", "17a", "18", "19", "20", "20a", "20b",
            "21", "22", "23", "24", "25",
            // Anticodon stem/loop (26-44)
            "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
            "39", "40", "41", "42", "43", "44",
            // Variable loop (45-47, variable length)
            "45", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e1", "e2", "e3", "e4", "e5",
            "e21", "e22", "e23", "e24", "e25", "e26", "e27", "46", "47",
            // T-stem/loop (48-65)
            "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61",
            "62", "63", "64", "65",
            // 3' acceptor stem
            "66", "67", "68", "69", "70", "71", "72",
            // Discriminator and CCA
            "73", "74", "75", "76",
        ];

        for (idx, pos) in positions.iter().enumerate() {
            let sprinzl = SprinzlPosition(pos.to_string());
            cm_to_sprinzl.insert(idx, sprinzl.clone());
            sprinzl_to_cm.insert(sprinzl, idx);
        }

        Self {
            cm_to_sprinzl,
            sprinzl_to_cm,
        }
    }

    /// Get Sprinzl position for a CM column index
    pub fn get_sprinzl(&self, cm_idx: usize) -> Option<&SprinzlPosition> {
        self.cm_to_sprinzl.get(&cm_idx)
    }

    /// Get CM column index for a Sprinzl position
    pub fn get_cm_index(&self, sprinzl: &SprinzlPosition) -> Option<usize> {
        self.sprinzl_to_cm.get(sprinzl).copied()
    }

    /// Map a sequence alignment to Sprinzl positions
    /// Returns a map from Sprinzl position to the sequence position (0-indexed)
    pub fn map_alignment(&self, alignment: &str) -> HashMap<SprinzlPosition, usize> {
        let mut result = HashMap::new();
        let mut seq_pos = 0;

        for (cm_idx, c) in alignment.chars().enumerate() {
            if c != '-' && c != '.' {
                if let Some(sprinzl) = self.cm_to_sprinzl.get(&cm_idx) {
                    result.insert(sprinzl.clone(), seq_pos);
                }
                seq_pos += 1;
            }
        }

        result
    }

    /// Check if a Sprinzl position is in a functionally important region
    pub fn is_critical_position(pos: &SprinzlPosition) -> bool {
        // Anticodon positions
        if ["34", "35", "36", "37"].contains(&pos.0.as_str()) {
            return true;
        }
        // D-loop modifications
        if ["16", "17", "20", "20a", "20b"].contains(&pos.0.as_str()) {
            return true;
        }
        // T-loop modifications
        if ["54", "55", "58"].contains(&pos.0.as_str()) {
            return true;
        }
        false
    }
}

impl Default for SprinzlMapper {
    fn default() -> Self {
        Self::new_standard()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_mapper() {
        let mapper = SprinzlMapper::new_standard();

        // Check some known positions
        assert!(mapper.get_sprinzl(0).is_some());
        assert_eq!(mapper.get_sprinzl(0).unwrap().0, "1");

        // Check anticodon wobble position (34)
        let pos34 = SprinzlPosition("34".to_string());
        assert!(mapper.get_cm_index(&pos34).is_some());
    }

    #[test]
    fn test_critical_positions() {
        assert!(SprinzlMapper::is_critical_position(&SprinzlPosition("34".to_string())));
        assert!(SprinzlMapper::is_critical_position(&SprinzlPosition("55".to_string())));
        assert!(!SprinzlMapper::is_critical_position(&SprinzlPosition("1".to_string())));
    }
}
