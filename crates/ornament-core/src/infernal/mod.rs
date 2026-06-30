//! Covariance-model interface module.
//!
//! Wraps the native `ornament-scfg` search engine. The `cmsearch` subprocess in
//! `runner` remains the differential-testing oracle and runtime fallback.

pub mod native;
pub mod parser;
pub mod runner;

pub use native::{scan_native, scan_native_multi};
pub use parser::{CMAlignment, CMHit};
pub use runner::InfernalRunner;
