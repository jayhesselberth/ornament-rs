//! `infernal-rs`: native Rust port of Infernal's covariance-model search internals.
//!
//! Modules (built incrementally):
//! - [`model`]  — the `CM`/states/nodes/transitions/emissions, sized from the alphabet's `K`.
//! - [`cmfile`] — the `.cm` v1.1 ASCII parser (the embedded HMMER3 HMM is retained as text).
//! - `config`   — `cm_Configure` equivalent (log-odds scores, local/glocal, consensus) [Phase 2/3].
//! - `search`   — non-banded CYK + Inside scanning DP [Phase 3].
//! - `align`    — CYK traceback -> parsetree -> aligned residues + RF/SS_cons [Phase 4].

pub mod align;
pub mod cmfile;
pub mod config;
pub mod emit;
pub mod emitmap;
pub mod evalues;
pub mod model;
pub mod search;

pub use align::{align_glocal, AlignedResidue, Alignment};
pub use cmfile::{parse_cm_file, parse_cm_str};
pub use config::{configure_local, configure_scores};
pub use emitmap::EmitMap;
pub use evalues::{evalue, score_to_evalue, SearchMode};
pub use model::{Cm, ExpInfo};
pub use search::{
    cyk_scan, cyk_scan_glocal, cyk_search, inside_scan, inside_scan_glocal, inside_search, CykHit,
    CykScan, Hit, Strand,
};

use thiserror::Error;

/// Errors produced by the Infernal CM layer.
#[derive(Debug, Error)]
pub enum InfernalError {
    #[error("I/O error: {0}")]
    Io(String),
    #[error("CM parse error: {0}")]
    Parse(String),
    #[error("easel error: {0}")]
    Easel(#[from] easel_rs::EaselError),
}
