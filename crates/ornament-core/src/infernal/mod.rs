//! Covariance-model interface module.
//!
//! Wraps the native `ornament-scfg` search engine. The `cmsearch` subprocess in
//! `runner` remains the differential-testing oracle and runtime fallback.

pub mod native;
pub mod parser;
pub mod report;
pub mod runner;
pub mod stockholm;

pub use native::{
    scan_native, scan_native_aligned, scan_native_multi, scan_native_multi_streaming,
};
pub use parser::{CMAlignment, CMHit};
pub use report::{write_tsv, write_tsv_rows, TSV_HEADER};
pub use runner::InfernalRunner;
pub use stockholm::{write_stockholm, AlignedRow, ModelMsa, ResCell};
