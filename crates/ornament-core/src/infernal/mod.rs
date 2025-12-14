//! Infernal interface module
//!
//! Provides wrappers around Infernal covariance model operations.

pub mod runner;
pub mod parser;

pub use runner::InfernalRunner;
pub use parser::{CMHit, CMAlignment};
