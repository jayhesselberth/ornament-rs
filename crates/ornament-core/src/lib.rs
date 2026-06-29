//! modtrna-core: Modification-aware tRNA analysis library
//!
//! This library provides tools for:
//! - Scanning genomic sequences for tRNAs using Infernal covariance models
//! - Mapping tRNA sequences to Sprinzl positions
//! - Analyzing modification compatibility at each position
//! - Detecting "odd" tRNAs with modification-incompatible variants

pub mod analysis;
pub mod infernal;
pub mod integration;
pub mod modification;
pub mod output;

// Re-export commonly used types
pub use analysis::TRNAHit;
pub use modification::types::{
    ConservationLevel, FunctionalRole, ModCode, Modification, PositionModExpectation, RnaBase,
    SprinzlPosition,
};
