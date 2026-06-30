//! E-value computation from the calibrated exponential-tail fits (`ExpInfo`).
//!
//! Faithful to Infernal's `Score2E`: a bit score is converted to an expected number of
//! chance hits using the per-mode exponential tail and the effective search-space size.
//! Verified against `cmsearch` output (the RF00005 consensus scores 87.3 bits → 1.1e-17
//! glocal CYK, 2e-18 glocal Inside).

use crate::model::{exp_mode, Cm, ExpInfo};

/// Search algorithm + model configuration, selecting which exponential tail to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    GlocalCyk,
    GlocalInside,
    LocalCyk,
    LocalInside,
}

impl SearchMode {
    fn exp_index(self) -> usize {
        match self {
            SearchMode::GlocalCyk => exp_mode::CM_GC,
            SearchMode::GlocalInside => exp_mode::CM_GI,
            SearchMode::LocalCyk => exp_mode::CM_LC,
            SearchMode::LocalInside => exp_mode::CM_LI,
        }
    }
}

/// E-value for a hit of bit score `score` over `searched_residues` total residues
/// (= sequence length for a single strand, `2*L` for both strands).
///
/// `E = (searched_residues / dbsize) * nrandhits * exp(-lambda * (score - mu_extrap))`.
pub fn score_to_evalue(exp: &ExpInfo, score: f32, searched_residues: f64) -> f64 {
    let cur_eff_dbsize = (searched_residues / exp.dbsize) * exp.nrandhits as f64;
    cur_eff_dbsize * (-exp.lambda * (score as f64 - exp.mu_extrap)).exp()
}

/// Convenience: E-value for a hit under the given search mode, using the CM's calibrated
/// tail for that mode.
pub fn evalue(cm: &Cm, mode: SearchMode, score: f32, searched_residues: f64) -> Option<f64> {
    let exp = &cm.exp[mode.exp_index()];
    if !exp.is_valid {
        return None;
    }
    Some(score_to_evalue(exp, score, searched_residues))
}
