//! Generic MSV filter (`p7_GMSV`) — the ungapped Multi-Segment Viterbi score, plus its
//! bit-score wrapper.
//!
//! Faithful port of HMMER's `generic_msv.c`. MSV is Viterbi restricted to ungapped local
//! single-segment alignments: there are no I/D states, the only core transition is
//! Mk→Mk+1 (with `tMM = 0`), B→Mk entry is uniform at `log(2/(M(M+1)))`, and an E→J/E→C
//! split governed by the expected-hit count `nu`. It carries its own length model
//! (`tloop`/`tmove` derived from L and `nu`) rather than the profile's N/C/J transitions —
//! only the match-emission scores are reused from the profile. Scores are natural-log
//! nats; the bit score is `(raw_nats - null_one(L)) / ln2`, matching `hmmsearch`.

use ornament_alphabet::{Alphabet, Dsq};

use crate::hmm::P7Hmm;
use crate::profile::{bg_freqs, null_one, P7Profile};

/// Matrix special-cell indices (HMMER `p7G_*`), as in [`crate::fwd`].
mod mx {
    pub const E: usize = 0;
    pub const N: usize = 1;
    pub const J: usize = 2;
    pub const B: usize = 3;
    pub const C: usize = 4;
}

const NINF: f32 = f32::NEG_INFINITY;

/// Expected number of hits (`nu` in `p7_GMSV`); 2.0 is HMMER's default for the MSV filter.
const NU: f32 = 2.0;

/// Max with `-inf` acting as the additive identity (the ⊕ of the max-plus semiring).
#[inline]
fn max(a: f32, b: f32) -> f32 {
    if a >= b {
        a
    } else {
        b
    }
}

/// Raw MSV score in nats for a 1-based digital sequence `dsq` (sentinels at 0 and L+1).
/// Uses `nu = 2.0`. The profile supplies only the match-emission scores; the B→Mk entry,
/// Mk→Mk+1, and length-model transitions are MSV's own.
pub fn msv_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    let m = prof.m;
    let l = dsq.len().saturating_sub(2);
    let stride = m + 1;
    let mut mmx = vec![NINF; (l + 1) * stride];
    let mut xmx = vec![NINF; (l + 1) * 5];

    // MSV-specific transitions (generic_msv.c).
    let tloop = (l as f32 / (l as f32 + 3.0)).ln(); // N/C/J loop
    let tmove = (3.0 / (l as f32 + 3.0)).ln(); // N/C/J move
    let tbmk = (2.0 / (m as f32 * (m as f32 + 1.0))).ln(); // uniform B->Mk entry
    let tej = ((NU - 1.0) / NU).ln(); // E->J
    let tec = (1.0 / NU).ln(); // E->C

    // Row 0.
    xmx[mx::N] = 0.0;
    xmx[mx::B] = tmove; // S->N->B, no N-tail
                        // E, C, J already NINF; M row 0 already NINF.

    for (i, &xb) in dsq.iter().enumerate().take(l + 1).skip(1) {
        let x = xb as usize;
        let msc = &prof.msc; // [Kp][M+1]
        let row = i * stride;
        let prow = (i - 1) * stride;
        let xrow = i * 5;
        let pxrow = (i - 1) * 5;

        mmx[row] = NINF;
        xmx[xrow + mx::E] = NINF;

        for k in 1..=m {
            // Ungapped: extend the diagonal (M_{k-1}, tMM=0) or start a new segment (B->Mk).
            mmx[row + k] = msc[x][k] + max(mmx[prow + k - 1], xmx[pxrow + mx::B] + tbmk);
            xmx[xrow + mx::E] = max(xmx[xrow + mx::E], mmx[row + k]);
        }

        // Specials (N/C/J emissions score 0; CC/NN/JJ are real loop scores here).
        xmx[xrow + mx::J] = max(xmx[pxrow + mx::J] + tloop, xmx[xrow + mx::E] + tej);
        xmx[xrow + mx::C] = max(xmx[pxrow + mx::C] + tloop, xmx[xrow + mx::E] + tec);
        xmx[xrow + mx::N] = xmx[pxrow + mx::N] + tloop;
        xmx[xrow + mx::B] = max(xmx[xrow + mx::N] + tmove, xmx[xrow + mx::J] + tmove);
    }

    xmx[l * 5 + mx::C] + tmove
}

/// MSV **bit score** for a digital sequence against the HMM (config + emissions + MSV DP +
/// null subtraction). MSV waives the core transition penalties and uses a flat entry, so —
/// unlike Viterbi vs Forward — it is *not* bounded by the other filters and is calibrated
/// on its own Gumbel scale (the `STATS LOCAL MSV` line).
pub fn msv_bits(hmm: &P7Hmm, abc: &Alphabet, dsq: &[Dsq]) -> f32 {
    let l = dsq.len().saturating_sub(2);
    let bg = bg_freqs(hmm.k);
    let prof = P7Profile::config_local(hmm, abc, &bg, l);
    let raw = msv_nats(&prof, dsq);
    (raw - null_one(l)) / std::f32::consts::LN_2
}
