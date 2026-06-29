//! Integration module
//!
//! Integration with external data sources like modkit.

pub mod modkit;

pub use modkit::{join_to_sprinzl, BedMethylRecord, ObservedModCall};
