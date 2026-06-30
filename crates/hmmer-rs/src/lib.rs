//! `hmmer-rs`: native Rust port of the HMMER p7 machinery Infernal uses as its CM
//! acceleration filter (`P7_HMM` / `P7_PROFILE` + MSV/Viterbi/Forward).
//!
//! The CM core (`infernal-rs`) is functionally complete without this — the p7 filter is an
//! acceleration stage that pre-selects windows for the expensive CYK/Inside DP. The final
//! reported scores still come from the CM, so the filter only affects speed.
//!
//! Modules:
//! - [`hmm`] — parser for the embedded HMMER3/f filter HMM (`fp7`).
//! - [`profile`] — the configured search profile (`P7_PROFILE`, multihit local).
//! - [`fwd`] — the shared full-model DP: Forward and Viterbi, one recursion over a semiring.
//! - [`vit`] — the Viterbi bit-score wrapper.
//! - [`msv`] — the ungapped MSV filter (`p7_GMSV`).
//! - [`pvalue`] — filter bit score → P-value via the calibrated Gumbel/exponential tails.
//!
//! `Viterbi ≤ Forward` always holds (same profile, max ≤ sum). MSV is a *different*, more
//! permissive model (it waives core transition penalties and uses a flat entry), so it is
//! not bounded by the other two and is calibrated on its own Gumbel scale.

pub mod fwd;
pub mod hmm;
pub mod logsum;
pub mod msv;
pub mod profile;
pub mod pvalue;
pub mod vit;

pub use fwd::{forward_bits, forward_nats, viterbi_nats};
pub use hmm::{parse_p7_hmm, EvParam, P7Hmm};
pub use msv::{msv_bits, msv_nats};
pub use profile::P7Profile;
pub use pvalue::{forward_pvalue, msv_pvalue, viterbi_pvalue};
pub use vit::viterbi_bits;

use thiserror::Error;

/// Errors from the HMMER p7 layer.
#[derive(Debug, Error)]
pub enum HmmerError {
    #[error("p7 parse error: {0}")]
    Parse(String),
}
