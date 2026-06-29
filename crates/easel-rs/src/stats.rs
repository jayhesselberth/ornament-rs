//! Extreme-value statistics used for E-value computation.
//!
//! Faithful ports of `esl_gumbel_surv` and `esl_exp_surv`. Infernal reports CM hit
//! significance from an exponential-tail fit (`ExpInfo`) and HMM-filter significance
//! from Gumbel fits; both reduce a bit score + fit parameters to a survival probability,
//! which is scaled by the search space to yield an E-value.

const SMALLX1: f64 = 5e-9;

/// Gumbel right-tail survival P(X > x) for location `mu`, scale `lambda`
/// (`esl_gumbel_surv`).
pub fn gumbel_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    let ey = -(-y).exp();
    // Use 1 - e^z ~ -z when e^-y is small, for numerical stability in the far tail.
    if ey.abs() < SMALLX1 {
        -ey
    } else {
        1.0 - ey.exp()
    }
}

/// Exponential right-tail survival P(X > x) for location `mu`, decay `lambda`
/// (`esl_exp_surv`). Returns 1.0 for `x < mu`.
pub fn exp_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    if x < mu {
        1.0
    } else {
        (-lambda * (x - mu)).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gumbel_surv_known_values() {
        // At x = mu, survival = 1 - e^{-1} ~ 0.6321.
        assert!((gumbel_surv(0.0, 0.0, 1.0) - (1.0 - (-1.0f64).exp())).abs() < 1e-12);
        // Far tail stays finite and monotone-decreasing.
        let a = gumbel_surv(10.0, 0.0, 1.0);
        let b = gumbel_surv(20.0, 0.0, 1.0);
        assert!(a > b && b > 0.0);
    }

    #[test]
    fn exp_surv_behaviour() {
        assert_eq!(exp_surv(-1.0, 0.0, 0.5), 1.0); // below mu
        assert!((exp_surv(0.0, 0.0, 0.5) - 1.0).abs() < 1e-12);
        assert!((exp_surv(2.0, 0.0, 0.5) - (-1.0f64).exp()).abs() < 1e-12);
    }
}
