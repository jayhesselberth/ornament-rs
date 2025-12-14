//! Core types for RNA modifications

use serde::{Deserialize, Serialize};
use std::fmt;

/// RNA nucleotide bases
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RnaBase {
    A,
    C,
    G,
    U,
}

impl RnaBase {
    /// Convert from DNA base character (T -> U)
    pub fn from_dna_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Self::A),
            'C' => Some(Self::C),
            'G' => Some(Self::G),
            'T' | 'U' => Some(Self::U),
            _ => None,
        }
    }

    /// Convert to RNA character
    pub fn to_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::U => 'U',
        }
    }

    /// Get the DNA equivalent character (U -> T)
    pub fn to_dna_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::U => 'T',
        }
    }

    /// Get the complement base
    pub fn complement(self) -> Self {
        match self {
            Self::A => Self::U,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::U => Self::A,
        }
    }
}

impl fmt::Display for RnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// Modification code representation
/// Supports multiple naming conventions from MODOMICS
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ModCode {
    /// Single character code (e.g., 'D' for dihydrouridine, 'Y' for pseudouridine)
    SingleChar(char),
    /// ChEBI ontology ID
    ChEBI(u32),
    /// MODOMICS unicode character (e.g., 'Ψ')
    Unicode(char),
    /// Short name (e.g., "m1A", "m5C", "t6A")
    ShortName(String),
}

impl fmt::Display for ModCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModCode::SingleChar(c) => write!(f, "{}", c),
            ModCode::ChEBI(id) => write!(f, "CHEBI:{}", id),
            ModCode::Unicode(c) => write!(f, "{}", c),
            ModCode::ShortName(s) => write!(f, "{}", s),
        }
    }
}

/// An RNA modification with its chemical and biological properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Modification {
    /// Full name (e.g., "pseudouridine", "1-methyladenosine")
    pub name: String,
    /// Short name (e.g., "Psi", "m1A")
    pub short_name: String,
    /// Primary code for identification
    pub code: ModCode,
    /// Alternative codes (for lookup)
    #[serde(default)]
    pub alt_codes: Vec<ModCode>,
    /// Parent nucleotide this modification is derived from
    pub parent_base: RnaBase,
    /// What we expect to see in genomic DNA/cDNA at this position
    /// (usually same as parent_base, but T instead of U)
    pub genomic_expectation: RnaBase,
    /// Bases that are incompatible with this modification
    /// (i.e., if we see these bases, this modification cannot exist)
    pub incompatible_bases: Vec<RnaBase>,
    /// ChEBI ID if known
    pub chebi_id: Option<u32>,
    /// MODOMICS unicode character if available
    pub modomics_unicode: Option<char>,
}

impl Modification {
    /// Check if an observed base is compatible with this modification
    pub fn is_compatible(&self, observed: RnaBase) -> bool {
        !self.incompatible_bases.contains(&observed)
    }

    /// Check if an observed base is the expected genomic base
    pub fn is_expected(&self, observed: RnaBase) -> bool {
        observed == self.genomic_expectation
    }
}

/// Conservation level of a modification across organisms
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConservationLevel {
    /// Found in virtually all tRNAs (e.g., Psi55, T54)
    Universal,
    /// Found in specific domains of life
    DomainSpecific,
    /// Found in specific isotypes only
    IsotypeSpecific,
    /// Uncommon or organism-specific
    Rare,
}

/// Functional role of a modification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FunctionalRole {
    /// Affects codon recognition (positions 34, 37)
    AnticodonFunction,
    /// Affects tRNA structure/stability
    StructuralStability,
    /// Identity element for aminoacyl-tRNA synthetase
    AminoacylationIdentity,
    /// Unknown or multiple roles
    Unknown,
}

/// Sprinzl tRNA position (1-76 with possible insertions like 17a, 20a, etc.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SprinzlPosition(pub String);

impl SprinzlPosition {
    pub fn new(pos: impl Into<String>) -> Self {
        Self(pos.into())
    }

    /// Create from numeric position
    pub fn from_num(n: u8) -> Self {
        Self(n.to_string())
    }

    /// Get the base position number (e.g., 17 from "17a")
    pub fn base_number(&self) -> Option<u8> {
        self.0
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect::<String>()
            .parse()
            .ok()
    }

    /// Check if this is an insertion position (e.g., "17a")
    pub fn is_insertion(&self) -> bool {
        self.0.chars().any(|c| c.is_ascii_alphabetic())
    }
}

impl fmt::Display for SprinzlPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<u8> for SprinzlPosition {
    fn from(n: u8) -> Self {
        Self::from_num(n)
    }
}

impl From<&str> for SprinzlPosition {
    fn from(s: &str) -> Self {
        Self::new(s)
    }
}

/// Expected modification at a specific tRNA position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionModExpectation {
    /// Sprinzl position
    pub position: SprinzlPosition,
    /// Expected modifications at this position (may be multiple alternatives)
    pub modifications: Vec<Modification>,
    /// How conserved this modification is
    pub conservation: ConservationLevel,
    /// Functional role
    pub functional_role: FunctionalRole,
    /// Which isotypes this applies to (empty = all)
    #[serde(default)]
    pub isotypes: Vec<String>,
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Self::Plus),
            '-' => Some(Self::Minus),
            _ => None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// tRNA isotype (amino acid specificity)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Isotype(pub String);

impl Isotype {
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }

    /// Standard isotypes
    pub const ALA: &'static str = "Ala";
    pub const ARG: &'static str = "Arg";
    pub const ASN: &'static str = "Asn";
    pub const ASP: &'static str = "Asp";
    pub const CYS: &'static str = "Cys";
    pub const GLN: &'static str = "Gln";
    pub const GLU: &'static str = "Glu";
    pub const GLY: &'static str = "Gly";
    pub const HIS: &'static str = "His";
    pub const ILE: &'static str = "Ile";
    pub const LEU: &'static str = "Leu";
    pub const LYS: &'static str = "Lys";
    pub const MET: &'static str = "Met";
    pub const PHE: &'static str = "Phe";
    pub const PRO: &'static str = "Pro";
    pub const SER: &'static str = "Ser";
    pub const THR: &'static str = "Thr";
    pub const TRP: &'static str = "Trp";
    pub const TYR: &'static str = "Tyr";
    pub const VAL: &'static str = "Val";
    pub const SEC: &'static str = "SeC"; // Selenocysteine
    pub const SUP: &'static str = "Sup"; // Suppressor
    pub const IMET: &'static str = "iMet"; // Initiator methionine
    pub const UNDET: &'static str = "Undet"; // Undetermined
}

impl fmt::Display for Isotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rna_base_from_dna() {
        assert_eq!(RnaBase::from_dna_char('A'), Some(RnaBase::A));
        assert_eq!(RnaBase::from_dna_char('T'), Some(RnaBase::U));
        assert_eq!(RnaBase::from_dna_char('U'), Some(RnaBase::U));
        assert_eq!(RnaBase::from_dna_char('N'), None);
    }

    #[test]
    fn test_sprinzl_position() {
        let pos = SprinzlPosition::from_num(34);
        assert_eq!(pos.base_number(), Some(34));
        assert!(!pos.is_insertion());

        let pos_ins = SprinzlPosition::new("17a");
        assert_eq!(pos_ins.base_number(), Some(17));
        assert!(pos_ins.is_insertion());
    }

    #[test]
    fn test_modification_compatibility() {
        let psi = Modification {
            name: "pseudouridine".to_string(),
            short_name: "Psi".to_string(),
            code: ModCode::Unicode('Ψ'),
            alt_codes: vec![ModCode::SingleChar('Y')],
            parent_base: RnaBase::U,
            genomic_expectation: RnaBase::U,
            incompatible_bases: vec![RnaBase::A, RnaBase::G, RnaBase::C],
            chebi_id: Some(17802),
            modomics_unicode: Some('Ψ'),
        };

        assert!(psi.is_compatible(RnaBase::U));
        assert!(!psi.is_compatible(RnaBase::A));
        assert!(!psi.is_compatible(RnaBase::G));
        assert!(psi.is_expected(RnaBase::U));
    }
}
