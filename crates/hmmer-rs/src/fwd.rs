//! Generic full-model DP (`p7_GForward` / `p7_GViterbi`) over a [`P7Profile`], plus the
//! bit-score wrappers.
//!
//! Faithful port of HMMER's `generic_fwdback.c` / `generic_viterbi.c`. The two DPs are
//! structurally identical — they differ only in how alternative sub-parses are combined
//! within a cell (the "⊕" operator): Forward sums them (log-sum-exp), Viterbi takes the
//! best (max). That single difference is captured by the [`Semiring`] parameter, mirroring
//! how `infernal-rs`'s scanning DP shares CYK/Inside. The "⊗" combiner stays addition in
//! both cases. Scores are natural-log nats; the final bit score is
//! `(raw_nats - null_one(L)) / ln2`, matching `hmmsearch`.

use easel_rs::{Alphabet, Dsq};

use crate::hmm::P7Hmm;
use crate::profile::{bg_freqs, null_one, p7p, xstate, P7Profile};

/// Matrix special-cell indices (HMMER `p7G_*`): note this differs from the profile's
/// `xsc` ordering (`p7P_*`, where C=3, B=4).
mod mx {
    pub const E: usize = 0;
    pub const N: usize = 1;
    pub const J: usize = 2;
    pub const B: usize = 3;
    pub const C: usize = 4;
}

const NINF: f32 = f32::NEG_INFINITY;

/// The "⊕" operation that combines alternative sub-parses within a DP cell. The identity
/// (empty disjunction) is `-inf` for both implementors.
pub(crate) trait Semiring {
    fn or(a: f32, b: f32) -> f32;
}

/// Log-sum-exp semiring → Forward (sum over all paths).
pub(crate) struct LogSum;
impl Semiring for LogSum {
    #[inline]
    fn or(a: f32, b: f32) -> f32 {
        flogsum(a, b)
    }
}

/// Max-plus semiring → Viterbi (best single path).
pub(crate) struct MaxPlus;
impl Semiring for MaxPlus {
    #[inline]
    fn or(a: f32, b: f32) -> f32 {
        if a >= b {
            a
        } else {
            b
        }
    }
}

/// Natural-log log-sum-exp: `log(e^a + e^b)`.
#[inline]
fn flogsum(a: f32, b: f32) -> f32 {
    if a == NINF {
        return b;
    }
    if b == NINF {
        return a;
    }
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    hi + ((lo - hi) as f64).exp().ln_1p() as f32
}

/// Raw Forward score in nats for a 1-based digital sequence `dsq` (with sentinels at 0 and
/// L+1). The profile must already be length-configured for L.
pub fn forward_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    dp_nats::<LogSum>(prof, dsq)
}

/// Raw Viterbi score in nats (best single path): the same matrix recursion as
/// [`forward_nats`] under the max-plus semiring (`p7_GViterbi`). [`crate::vit`] has the
/// bit-score wrapper.
pub fn viterbi_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    dp_nats::<MaxPlus>(prof, dsq)
}

/// Shared Forward/Viterbi recursion, parameterized by the within-cell ⊕ ([`Semiring`]).
/// With `S = LogSum` this is `p7_GForward`; with `S = MaxPlus`, `p7_GViterbi`. The Dk→E
/// exits are folded in at every k; that is required for Forward and harmless for Viterbi
/// (a Dk→E path can never beat the best Mj→E in local/glocal mode, as noted in
/// `generic_viterbi.c`).
fn dp_nats<S: Semiring>(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    let m = prof.m;
    let l = dsq.len().saturating_sub(2);
    let stride = m + 1;
    let mut mmx = vec![NINF; (l + 1) * stride];
    let mut imx = vec![NINF; (l + 1) * stride];
    let mut dmx = vec![NINF; (l + 1) * stride];
    let mut xmx = vec![NINF; (l + 1) * 5];

    let esc = 0.0f32; // local mode: E exit from any match/delete
    let tsc = |ty: usize, k: usize| prof.tsc(ty, k);

    // Row 0.
    xmx[mx::N] = 0.0;
    xmx[mx::B] = prof.xsc[xstate::N][xstate::MOVE];
    // E, C, J already NINF; M/I/D row 0 already NINF.

    for (i, &xb) in dsq.iter().enumerate().take(l + 1).skip(1) {
        let x = xb as usize;
        let msc = &prof.msc; // [Kp][M+1]
        let isc = &prof.isc;
        let row = i * stride;
        let prow = (i - 1) * stride;
        let xrow = i * 5;
        let pxrow = (i - 1) * 5;

        mmx[row] = NINF;
        imx[row] = NINF;
        dmx[row] = NINF;
        xmx[xrow + mx::E] = NINF;

        for k in 1..m {
            // match
            let sc = S::or(
                S::or(
                    mmx[prow + k - 1] + tsc(p7p::MM, k - 1),
                    imx[prow + k - 1] + tsc(p7p::IM, k - 1),
                ),
                S::or(
                    xmx[pxrow + mx::B] + tsc(p7p::BM, k - 1),
                    dmx[prow + k - 1] + tsc(p7p::DM, k - 1),
                ),
            );
            mmx[row + k] = sc + msc[x][k];

            // insert
            let sc = S::or(
                mmx[prow + k] + tsc(p7p::MI, k),
                imx[prow + k] + tsc(p7p::II, k),
            );
            imx[row + k] = sc + isc[x][k];

            // delete
            dmx[row + k] = S::or(
                mmx[row + k - 1] + tsc(p7p::MD, k - 1),
                dmx[row + k - 1] + tsc(p7p::DD, k - 1),
            );

            // E
            xmx[xrow + mx::E] = S::or(
                S::or(mmx[row + k] + esc, dmx[row + k] + esc),
                xmx[xrow + mx::E],
            );
        }

        // unrolled M_M
        let sc = S::or(
            S::or(
                mmx[prow + m - 1] + tsc(p7p::MM, m - 1),
                imx[prow + m - 1] + tsc(p7p::IM, m - 1),
            ),
            S::or(
                xmx[pxrow + mx::B] + tsc(p7p::BM, m - 1),
                dmx[prow + m - 1] + tsc(p7p::DM, m - 1),
            ),
        );
        mmx[row + m] = sc + msc[x][m];
        imx[row + m] = NINF;
        dmx[row + m] = S::or(
            mmx[row + m - 1] + tsc(p7p::MD, m - 1),
            dmx[row + m - 1] + tsc(p7p::DD, m - 1),
        );
        xmx[xrow + mx::E] = S::or(S::or(mmx[row + m], dmx[row + m]), xmx[xrow + mx::E]);

        // specials
        xmx[xrow + mx::J] = S::or(
            xmx[pxrow + mx::J] + prof.xsc[xstate::J][xstate::LOOP],
            xmx[xrow + mx::E] + prof.xsc[xstate::E][xstate::LOOP],
        );
        xmx[xrow + mx::C] = S::or(
            xmx[pxrow + mx::C] + prof.xsc[xstate::C][xstate::LOOP],
            xmx[xrow + mx::E] + prof.xsc[xstate::E][xstate::MOVE],
        );
        xmx[xrow + mx::N] = xmx[pxrow + mx::N] + prof.xsc[xstate::N][xstate::LOOP];
        xmx[xrow + mx::B] = S::or(
            xmx[xrow + mx::N] + prof.xsc[xstate::N][xstate::MOVE],
            xmx[xrow + mx::J] + prof.xsc[xstate::J][xstate::MOVE],
        );
    }

    xmx[l * 5 + mx::C] + prof.xsc[xstate::C][xstate::MOVE]
}

/// Forward **bit score** for a digital sequence against the HMM (config + length model +
/// Forward + null subtraction), matching `hmmsearch`'s full-sequence Forward score.
pub fn forward_bits(hmm: &P7Hmm, abc: &Alphabet, dsq: &[Dsq]) -> f32 {
    let l = dsq.len().saturating_sub(2);
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(hmm, abc, &bg, l);
    let raw = forward_nats(&prof, dsq);
    (raw - null_one(l)) / std::f32::consts::LN_2
}
