//! Modification database - MODOMICS-derived tRNA modification expectations

use crate::modification::types::*;
use rustc_hash::FxHashMap;

/// Database of known tRNA modifications and their position expectations
pub struct ModificationDatabase {
    /// All known modifications indexed by short name
    modifications: FxHashMap<String, Modification>,
    /// Position-specific expectations (Sprinzl position -> expectations)
    position_expectations: FxHashMap<String, Vec<PositionModExpectation>>,
}

impl ModificationDatabase {
    /// Create a new database with default eukaryotic modifications
    pub fn eukaryotic() -> Self {
        let mut db = Self {
            modifications: FxHashMap::default(),
            position_expectations: FxHashMap::default(),
        };
        db.load_default_modifications();
        db.load_eukaryotic_expectations();
        db
    }

    /// Get a modification by short name
    pub fn get_modification(&self, short_name: &str) -> Option<&Modification> {
        self.modifications.get(short_name)
    }

    /// Get all modifications in the database
    pub fn modifications(&self) -> &FxHashMap<String, Modification> {
        &self.modifications
    }

    /// Get all expected modifications at a Sprinzl position
    pub fn get_expectations(&self, position: &SprinzlPosition) -> Vec<&PositionModExpectation> {
        self.position_expectations
            .get(&position.0)
            .map(|v| v.iter().collect())
            .unwrap_or_default()
    }

    /// Get expectations for a specific isotype at a position
    pub fn get_expectations_for_isotype(
        &self,
        position: &SprinzlPosition,
        isotype: &Isotype,
    ) -> Vec<&PositionModExpectation> {
        self.get_expectations(position)
            .into_iter()
            .filter(|exp| exp.isotypes.is_empty() || exp.isotypes.contains(&isotype.0))
            .collect()
    }

    fn load_default_modifications(&mut self) {
        // Pseudouridine (Psi/Y) - most common modification
        self.add_modification(Modification {
            name: "pseudouridine".to_string(),
            short_name: "Psi".to_string(),
            code: ModCode::Unicode('Ψ'),
            alt_codes: vec![ModCode::SingleChar('Y')],
            parent_base: RnaBase::U,
            genomic_expectation: RnaBase::U,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::C],
            chebi_id: Some(17802),
            modomics_unicode: Some('Ψ'),
        });

        // Dihydrouridine (D)
        self.add_modification(Modification {
            name: "dihydrouridine".to_string(),
            short_name: "D".to_string(),
            code: ModCode::SingleChar('D'),
            alt_codes: vec![],
            parent_base: RnaBase::U,
            genomic_expectation: RnaBase::U,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::C],
            chebi_id: Some(15802),
            modomics_unicode: Some('D'),
        });

        // 5-methyluridine (m5U/T/rT) - ribothymidine
        self.add_modification(Modification {
            name: "5-methyluridine".to_string(),
            short_name: "m5U".to_string(),
            code: ModCode::ShortName("m5U".to_string()),
            alt_codes: vec![ModCode::SingleChar('T')],
            parent_base: RnaBase::U,
            genomic_expectation: RnaBase::U,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::C],
            chebi_id: Some(16695),
            modomics_unicode: Some('T'),
        });

        // 1-methyladenosine (m1A)
        self.add_modification(Modification {
            name: "1-methyladenosine".to_string(),
            short_name: "m1A".to_string(),
            code: ModCode::ShortName("m1A".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::A,
            genomic_expectation: RnaBase::A,
            incompatible_bases: vec![RnaBase::G, RnaBase::C, RnaBase::U],
            chebi_id: Some(21837),
            modomics_unicode: Some('"'),
        });

        // 1-methylguanosine (m1G)
        self.add_modification(Modification {
            name: "1-methylguanosine".to_string(),
            short_name: "m1G".to_string(),
            code: ModCode::ShortName("m1G".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::G,
            genomic_expectation: RnaBase::G,
            incompatible_bases: vec![RnaBase::A, RnaBase::C, RnaBase::U],
            chebi_id: Some(21836),
            modomics_unicode: Some('K'),
        });

        // N6-threonylcarbamoyladenosine (t6A)
        self.add_modification(Modification {
            name: "N6-threonylcarbamoyladenosine".to_string(),
            short_name: "t6A".to_string(),
            code: ModCode::ShortName("t6A".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::A,
            genomic_expectation: RnaBase::A,
            incompatible_bases: vec![RnaBase::G, RnaBase::C, RnaBase::U],
            chebi_id: Some(20817),
            modomics_unicode: Some('6'),
        });

        // N6-isopentenyladenosine (i6A)
        self.add_modification(Modification {
            name: "N6-isopentenyladenosine".to_string(),
            short_name: "i6A".to_string(),
            code: ModCode::ShortName("i6A".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::A,
            genomic_expectation: RnaBase::A,
            incompatible_bases: vec![RnaBase::G, RnaBase::C, RnaBase::U],
            chebi_id: Some(17588),
            modomics_unicode: Some('+'),
        });

        // Inosine (I) - A to I editing at wobble position
        self.add_modification(Modification {
            name: "inosine".to_string(),
            short_name: "I".to_string(),
            code: ModCode::SingleChar('I'),
            alt_codes: vec![],
            parent_base: RnaBase::A,
            genomic_expectation: RnaBase::A,
            incompatible_bases: vec![RnaBase::G, RnaBase::C, RnaBase::U],
            chebi_id: Some(17596),
            modomics_unicode: Some('I'),
        });

        // Queuosine (Q)
        self.add_modification(Modification {
            name: "queuosine".to_string(),
            short_name: "Q".to_string(),
            code: ModCode::SingleChar('Q'),
            alt_codes: vec![],
            parent_base: RnaBase::G,
            genomic_expectation: RnaBase::G,
            incompatible_bases: vec![RnaBase::A, RnaBase::C, RnaBase::U],
            chebi_id: Some(17399),
            modomics_unicode: Some('Q'),
        });

        // 2'-O-methylcytidine (Cm)
        self.add_modification(Modification {
            name: "2'-O-methylcytidine".to_string(),
            short_name: "Cm".to_string(),
            code: ModCode::ShortName("Cm".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::C,
            genomic_expectation: RnaBase::C,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::U],
            chebi_id: Some(19228),
            modomics_unicode: Some('B'),
        });

        // 5-methylcytidine (m5C)
        self.add_modification(Modification {
            name: "5-methylcytidine".to_string(),
            short_name: "m5C".to_string(),
            code: ModCode::ShortName("m5C".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::C,
            genomic_expectation: RnaBase::C,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::U],
            chebi_id: Some(27480),
            modomics_unicode: Some('?'),
        });

        // 7-methylguanosine (m7G)
        self.add_modification(Modification {
            name: "7-methylguanosine".to_string(),
            short_name: "m7G".to_string(),
            code: ModCode::ShortName("m7G".to_string()),
            alt_codes: vec![],
            parent_base: RnaBase::G,
            genomic_expectation: RnaBase::G,
            incompatible_bases: vec![RnaBase::A, RnaBase::C, RnaBase::U],
            chebi_id: Some(2274),
            modomics_unicode: Some('7'),
        });
    }

    fn load_eukaryotic_expectations(&mut self) {
        // Position 8 - Usually A, sometimes modified
        // (Not adding modification requirement)

        // Position 13-17 - D-loop dihydrouridines
        for pos in [16, 17, 20] {
            self.add_position_expectation(PositionModExpectation {
                position: SprinzlPosition::from_num(pos),
                modifications: vec![self.modifications.get("D").unwrap().clone()],
                conservation: ConservationLevel::Universal,
                functional_role: FunctionalRole::StructuralStability,
                isotypes: vec![],
            });
        }

        // Position 32 - Cm in some tRNAs
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(32),
            modifications: vec![self.modifications.get("Cm").unwrap().clone()],
            conservation: ConservationLevel::IsotypeSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![Isotype::PHE.to_string(), Isotype::TRP.to_string()],
        });

        // Position 34 - Wobble position - various modifications
        // Inosine in specific tRNAs
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(34),
            modifications: vec![self.modifications.get("I").unwrap().clone()],
            conservation: ConservationLevel::IsotypeSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![
                Isotype::ALA.to_string(),
                Isotype::ARG.to_string(),
                Isotype::ILE.to_string(),
                Isotype::LEU.to_string(),
                Isotype::PRO.to_string(),
                Isotype::SER.to_string(),
                Isotype::THR.to_string(),
                Isotype::VAL.to_string(),
            ],
        });

        // Queuosine at position 34 for specific isotypes
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(34),
            modifications: vec![self.modifications.get("Q").unwrap().clone()],
            conservation: ConservationLevel::IsotypeSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![
                Isotype::ASN.to_string(),
                Isotype::ASP.to_string(),
                Isotype::HIS.to_string(),
                Isotype::TYR.to_string(),
            ],
        });

        // Position 37 - 3' of anticodon - hypermodified in most tRNAs
        // t6A is common
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(37),
            modifications: vec![self.modifications.get("t6A").unwrap().clone()],
            conservation: ConservationLevel::DomainSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![
                Isotype::ILE.to_string(),
                Isotype::LYS.to_string(),
                Isotype::ASN.to_string(),
                Isotype::SER.to_string(),
                Isotype::THR.to_string(),
            ],
        });

        // i6A at position 37
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(37),
            modifications: vec![self.modifications.get("i6A").unwrap().clone()],
            conservation: ConservationLevel::IsotypeSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![
                Isotype::CYS.to_string(),
                Isotype::SER.to_string(),
                Isotype::TRP.to_string(),
            ],
        });

        // m1G at position 37 for some isotypes
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(37),
            modifications: vec![self.modifications.get("m1G").unwrap().clone()],
            conservation: ConservationLevel::IsotypeSpecific,
            functional_role: FunctionalRole::AnticodonFunction,
            isotypes: vec![
                Isotype::ALA.to_string(),
                Isotype::ARG.to_string(),
                Isotype::LEU.to_string(),
                Isotype::PRO.to_string(),
            ],
        });

        // Position 46 - m7G
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(46),
            modifications: vec![self.modifications.get("m7G").unwrap().clone()],
            conservation: ConservationLevel::Universal,
            functional_role: FunctionalRole::StructuralStability,
            isotypes: vec![],
        });

        // Position 48/49 - m5C
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(48),
            modifications: vec![self.modifications.get("m5C").unwrap().clone()],
            conservation: ConservationLevel::DomainSpecific,
            functional_role: FunctionalRole::StructuralStability,
            isotypes: vec![],
        });

        // Position 54 - m5U (ribothymidine) - nearly universal
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(54),
            modifications: vec![self.modifications.get("m5U").unwrap().clone()],
            conservation: ConservationLevel::Universal,
            functional_role: FunctionalRole::StructuralStability,
            isotypes: vec![],
        });

        // Position 55 - Pseudouridine - universal
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(55),
            modifications: vec![self.modifications.get("Psi").unwrap().clone()],
            conservation: ConservationLevel::Universal,
            functional_role: FunctionalRole::StructuralStability,
            isotypes: vec![],
        });

        // Position 58 - m1A - very common
        self.add_position_expectation(PositionModExpectation {
            position: SprinzlPosition::from_num(58),
            modifications: vec![self.modifications.get("m1A").unwrap().clone()],
            conservation: ConservationLevel::Universal,
            functional_role: FunctionalRole::StructuralStability,
            isotypes: vec![],
        });
    }

    fn add_modification(&mut self, modification: Modification) {
        self.modifications
            .insert(modification.short_name.clone(), modification);
    }

    fn add_position_expectation(&mut self, expectation: PositionModExpectation) {
        self.position_expectations
            .entry(expectation.position.0.clone())
            .or_default()
            .push(expectation);
    }
}

impl Default for ModificationDatabase {
    fn default() -> Self {
        Self::eukaryotic()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_database_creation() {
        let db = ModificationDatabase::eukaryotic();

        // Should have pseudouridine
        assert!(db.get_modification("Psi").is_some());

        // Position 55 should have Psi expectation
        let exp55 = db.get_expectations(&SprinzlPosition::from_num(55));
        assert!(!exp55.is_empty());
        assert_eq!(exp55[0].modifications[0].short_name, "Psi");
    }

    #[test]
    fn test_isotype_specific_expectations() {
        let db = ModificationDatabase::eukaryotic();

        // Position 34 for Ala should have Inosine
        let exp34_ala = db.get_expectations_for_isotype(
            &SprinzlPosition::from_num(34),
            &Isotype::new("Ala"),
        );
        assert!(exp34_ala.iter().any(|e|
            e.modifications.iter().any(|m| m.short_name == "I")
        ));
    }
}
