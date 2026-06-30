//! Viterbi filter bit-score wrapper.
//!
//! The Viterbi DP itself lives in [`crate::fwd`] — it is the shared full-model recursion
//! under the max-plus semiring (`p7_GViterbi`). Here we only add the profile config +
//! length model + null subtraction, mirroring [`crate::fwd::forward_bits`].

use ornament_alphabet::{Alphabet, Dsq};

use crate::fwd::viterbi_nats;
use crate::hmm::P7Hmm;
use crate::profile::{bg_freqs, null_one, P7Profile};

/// Viterbi **bit score** for a digital sequence against the HMM (config + length model +
/// Viterbi + null subtraction). Uses the same multihit-local profile as Forward, so the
/// optimal-path score is always ≤ the Forward (sum-over-paths) score.
pub fn viterbi_bits(hmm: &P7Hmm, abc: &Alphabet, dsq: &[Dsq]) -> f32 {
    let l = dsq.len().saturating_sub(2);
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(hmm, abc, &bg, l);
    let raw = viterbi_nats(&prof, dsq);
    (raw - null_one(l)) / std::f32::consts::LN_2
}
