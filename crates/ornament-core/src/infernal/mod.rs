//! Infernal interface module
//!
//! Provides wrappers around Infernal covariance model operations.

pub mod ffi;
pub mod runner;
pub mod parser;

pub use ffi::{Alphabet, CovarianceModel, HmmFilter, Sequence, SequenceFile, TopHits};
pub use runner::InfernalRunner;
pub use parser::{CMHit, CMAlignment};
