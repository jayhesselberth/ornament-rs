//! Filter-score P-values from the calibrated extreme-value fits (`STATS LOCAL` lines).
//!
//! Faithful to HMMER's pipeline (`p7_pli.c`): a filter **bit** score is converted to a
//! P-value through that filter's calibrated right-tail —
//!
//! - **MSV** and **Viterbi** use a **Gumbel** survival (`esl_gumbel_surv`), and
//! - **Forward** uses an **exponential** tail (`esl_exp_surv`),
//!
//! exactly as `p7_Pipeline` does (`P = esl_gumbel_surv(seq_score, MMU, MLAMBDA)`,
//! `P = esl_exp_surv(seq_score, FTAU, FLAMBDA)`).
//!
//! Units: the `STATS LOCAL` `mu`/`lambda` are calibrated against the **bit** score (the
//! same `(raw_nats - null_one(L)) / ln2` that [`crate::forward_bits`] etc. return — note
//! `lambda ~ 0.69 = ln2`, i.e. already bit-space). The bit score's length normalization
//! (length-model reconfig + null subtraction) makes its distribution approximately
//! length-independent, so HMMER applies the single stored `mu`/`lambda` directly, with no
//! per-length correction — and so do we.

use easel_rs::stats::{exp_surv, gumbel_surv};

use crate::hmm::P7Hmm;

/// MSV-filter P-value for an MSV **bit** score: Gumbel right-tail under the model's
/// `STATS LOCAL MSV mu lambda`.
pub fn msv_pvalue(hmm: &P7Hmm, bits: f32) -> f64 {
    gumbel_surv(bits as f64, hmm.msv.mu, hmm.msv.lambda)
}

/// Viterbi-filter P-value for a Viterbi **bit** score: Gumbel right-tail under the model's
/// `STATS LOCAL VITERBI mu lambda`.
pub fn viterbi_pvalue(hmm: &P7Hmm, bits: f32) -> f64 {
    gumbel_surv(bits as f64, hmm.vit.mu, hmm.vit.lambda)
}

/// Forward-filter P-value for a Forward **bit** score: exponential right-tail under the
/// model's `STATS LOCAL FORWARD tau lambda` (HMMER fits Forward's tail to an exponential,
/// not a Gumbel).
pub fn forward_pvalue(hmm: &P7Hmm, bits: f32) -> f64 {
    exp_surv(bits as f64, hmm.fwd.mu, hmm.fwd.lambda)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmm::EvParam;

    /// A stand-in HMM carrying only the calibration params (the RF00005 tRNA `STATS` lines).
    fn trna_stats() -> P7Hmm {
        P7Hmm {
            name: "tRNA".into(),
            acc: None,
            m: 71,
            k: 4,
            maxl: 253,
            mat: Vec::new(),
            ins: Vec::new(),
            t: Vec::new(),
            compo: None,
            msv: EvParam {
                mu: -8.9378,
                lambda: 0.73469,
            },
            vit: EvParam {
                mu: -9.6821,
                lambda: 0.73469,
            },
            fwd: EvParam {
                mu: -2.5854,
                lambda: 0.73469,
            },
        }
    }

    /// At `bits = mu`, the Gumbel survival is exactly `1 - e^{-1} ≈ 0.6321` — the hand
    /// anchor for "a score right at the location parameter is a coin-flip-ish P".
    #[test]
    fn gumbel_pvalue_at_mu_is_063() {
        let hmm = trna_stats();
        let want = 1.0 - (-1.0f64).exp();
        assert!((msv_pvalue(&hmm, hmm.msv.mu as f32) - want).abs() < 1e-4);
        assert!((viterbi_pvalue(&hmm, hmm.vit.mu as f32) - want).abs() < 1e-4);
    }

    /// A strong filter score must yield a tiny P-value. The validated 53.46-bit Forward hit
    /// (vs hmmsearch's 53.5) gives `exp(-0.73469 * (53.46 + 2.5854)) ≈ 1.3e-18`.
    #[test]
    fn forward_pvalue_high_score_is_tiny() {
        let hmm = trna_stats();
        let p = forward_pvalue(&hmm, 53.46);
        // Cross-check against the closed form.
        let want = (-0.73469f64 * (53.46 + 2.5854)).exp();
        assert!(
            (p / want - 1.0).abs() < 1e-6,
            "p = {p:.3e}, want {want:.3e}"
        );
        assert!(p < 1e-15, "53-bit Forward P should be tiny, got {p:.3e}");
        // The exponential tail is flat below tau (returns 1.0 for x < tau).
        assert_eq!(forward_pvalue(&hmm, -10.0), 1.0);
    }

    /// MSV/Viterbi P-values are monotone-decreasing in the score and small in the tail.
    #[test]
    fn gumbel_pvalues_monotone_and_small() {
        let hmm = trna_stats();
        let hi = msv_pvalue(&hmm, 40.0);
        let lo = msv_pvalue(&hmm, 50.0);
        assert!(hi > lo && lo > 0.0, "monotone: {hi:.3e} > {lo:.3e}");
        assert!(msv_pvalue(&hmm, 40.0) < 1e-12);
        assert!(viterbi_pvalue(&hmm, 40.0) < 1e-12);
    }
}
