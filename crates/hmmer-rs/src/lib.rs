//! `hmmer-rs`: native Rust port of the HMMER p7 machinery Infernal uses as its CM
//! acceleration filter (`P7_HMM` / `P7_PROFILE` + MSV/Viterbi/Forward).
//!
//! The CM core (`infernal-rs`) is functionally complete without this — the p7 filter is an
//! acceleration stage that pre-selects windows for the expensive CYK/Inside DP. The final
//! reported scores still come from the CM, so the filter only affects speed.
//!
//! Modules:
//! - [`hmm`] — parser for the embedded HMMER3/f filter HMM (`fp7`).
//! - `profile`, `msv`, `vit`, `fwd` — to come.

pub mod fwd;
pub mod hmm;
pub mod profile;

pub use fwd::{forward_bits, forward_nats};
pub use hmm::{parse_p7_hmm, EvParam, P7Hmm};
pub use profile::P7Profile;

use thiserror::Error;

/// Errors from the HMMER p7 layer.
#[derive(Debug, Error)]
pub enum HmmerError {
    #[error("p7 parse error: {0}")]
    Parse(String),
}
