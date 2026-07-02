//! `ornament-alphabet`: the modification-aware digital alphabet and digital
//! sequences underpinning Ornament's covariance-model search.
//!
//! Originally a native Rust port of the [Easel](https://github.com/EddyRivasLab/easel)
//! alphabet subset, **extended** with first-class modified-base symbols. FASTA input is
//! provided by the [`noodles`](https://github.com/zaeleus/noodles) library.
//!
//! Modules:
//! - [`alphabet`] — digital alphabets with first-class modified-base symbols.
//! - [`sq`] — sequences and FASTA input.
//!
//! (Extreme-value E-value statistics live in the separate `ornament-stats` crate.)

pub mod alphabet;
pub mod sq;

pub use alphabet::{Alphabet, AlphabetDescriptor, AlphabetType, Dsq, ModSpec, ModSymbol};
pub use sq::{read_fasta, read_fasta_from, Sequence};

use thiserror::Error;

/// Errors produced by the alphabet / sequence layer.
#[derive(Debug, Error)]
pub enum AlphabetError {
    #[error("alphabet error: {0}")]
    Alphabet(String),
    #[error("I/O error: {0}")]
    Io(String),
    #[error("parse error: {0}")]
    Parse(String),
}
