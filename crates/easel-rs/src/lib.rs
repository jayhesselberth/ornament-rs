//! `easel-rs`: a native Rust port of the [Easel](https://github.com/EddyRivasLab/easel)
//! subset used by Infernal's covariance-model search, **extended** with a
//! modification-aware alphabet.
//!
//! Modules:
//! - [`alphabet`] — digital alphabets with first-class modified-base symbols.
//! - [`sq`] — sequences and FASTA input.
//! - [`stats`] — Gumbel / exponential survival functions for E-values.

pub mod alphabet;
pub mod sq;
pub mod stats;

pub use alphabet::{Alphabet, AlphabetType, Dsq, ModSpec, ModSymbol};
pub use sq::{read_fasta, read_fasta_from, Sequence};

use thiserror::Error;

/// Errors produced by the Easel layer.
#[derive(Debug, Error)]
pub enum EaselError {
    #[error("alphabet error: {0}")]
    Alphabet(String),
    #[error("I/O error: {0}")]
    Io(String),
    #[error("parse error: {0}")]
    Parse(String),
}
