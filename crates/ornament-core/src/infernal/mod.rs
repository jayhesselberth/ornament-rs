//! Infernal interface module
//!
//! Provides wrappers around Infernal covariance model operations.

// Legacy C FFI wrappers — only compiled with the `ffi` feature. The native Rust port
// (easel-rs / infernal-rs) replaces these; the `cmsearch` subprocess in `runner` remains
// the differential-testing oracle and runtime fallback.
#[cfg(feature = "ffi")]
pub mod ffi;
pub mod runner;
pub mod parser;

#[cfg(feature = "ffi")]
pub use ffi::{Alphabet, CovarianceModel, HmmFilter, Sequence, SequenceFile, TopHits};
pub use runner::InfernalRunner;
pub use parser::{CMHit, CMAlignment};
