//! `hmmer-rs`: native Rust port of the HMMER p7 machinery Infernal uses as its CM
//! acceleration filter (`P7_HMM`/`P7_PROFILE`/`P7_OPROFILE`/`P7_BG` + MSV/Viterbi/Forward).
//!
//! Scaffolding only — implemented in Phase 5. The correctness-first CM core (Phases 1–4)
//! runs the unfiltered CYK/Inside DP and does not depend on this crate yet.

/// Placeholder until the p7 filter pipeline lands in Phase 5.
pub fn placeholder() {}
