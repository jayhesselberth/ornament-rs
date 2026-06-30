//! Probability-space (odds-ratio) Forward filter — the representation HMMER's striped-SIMD
//! `p7_ForwardFilter` uses (`impl_sse/fwdback.c`). Unlike the log-space [`crate::fwd`] DP
//! (which needs a log-sum-exp per cell), the recursion here is pure multiply-add on odds
//! ratios, so it vectorizes cleanly. Underflow is held off by **sparse row rescaling**: when
//! the row's E odds exceed a threshold the whole row is divided down and `log(scale)` is
//! accumulated, exactly as HMMER does. The returned score is back in nats:
//! `totscale + log(xC · t_C,move)`, matching `forward_nats` to filter precision.
//!
//! This module is built in stages: this file is the scalar (non-striped) odds-space engine
//! that establishes correctness and serves as the oracle for the striped SIMD kernel.

use easel_rs::Dsq;

use crate::profile::{p7p, xstate, P7Profile};

/// Sparse-rescaling trigger (HMMER's `forward_engine`): once the row's summed E odds exceed
/// this (~e^10, ~10% of f32 dynamic range) the row is divided by `xE` and `log(xE)` banked.
const RESCALE_THRESHOLD: f32 = 1.0e4;

/// Precomputed odds ratios (`exp` of the log-odds profile), the multiply-add form the Forward
/// filter consumes. Transition arrays are indexed by node `k` exactly as [`P7Profile::tsc`].
struct OddsProfile {
    /// Match-emission odds per residue code, `me[x][k] = exp(msc[x][k])`.
    me: Vec<Vec<f32>>,
    // Transition odds, indexed by the same `k` as `prof.tsc(_, k)`.
    tbm: Vec<f32>,
    tmm: Vec<f32>,
    tim: Vec<f32>,
    tdm: Vec<f32>,
    tmd: Vec<f32>,
    tdd: Vec<f32>,
    tmi: Vec<f32>,
    tii: Vec<f32>,
    // Special-state odds.
    n_loop: f32,
    n_move: f32,
    c_loop: f32,
    c_move: f32,
    j_loop: f32,
    j_move: f32,
    e_move: f32,
    e_loop: f32,
}

impl OddsProfile {
    fn build(prof: &P7Profile) -> OddsProfile {
        let m = prof.m;
        let exp_t = |ty: usize| (0..=m).map(|k| prof.tsc(ty, k).exp()).collect::<Vec<f32>>();
        let me = prof
            .msc
            .iter()
            .map(|row| row.iter().map(|&s| s.exp()).collect())
            .collect();
        let xf = |st: usize, dir: usize| prof.xsc[st][dir].exp();
        OddsProfile {
            me,
            tbm: exp_t(p7p::BM),
            tmm: exp_t(p7p::MM),
            tim: exp_t(p7p::IM),
            tdm: exp_t(p7p::DM),
            tmd: exp_t(p7p::MD),
            tdd: exp_t(p7p::DD),
            tmi: exp_t(p7p::MI),
            tii: exp_t(p7p::II),
            n_loop: xf(xstate::N, xstate::LOOP),
            n_move: xf(xstate::N, xstate::MOVE),
            c_loop: xf(xstate::C, xstate::LOOP),
            c_move: xf(xstate::C, xstate::MOVE),
            j_loop: xf(xstate::J, xstate::LOOP),
            j_move: xf(xstate::J, xstate::MOVE),
            e_move: xf(xstate::E, xstate::MOVE),
            e_loop: xf(xstate::E, xstate::LOOP),
        }
    }
}

/// Raw Forward score in nats via the odds-space DP with sparse rescaling. Matches
/// [`crate::fwd::forward_nats`] to filter precision (the difference is f32 odds-space rounding,
/// well under the filter's bit thresholds). The profile must be length-configured for `L`.
#[allow(clippy::needless_range_loop)] // DP indices (k, i) address several parallel rows
pub fn forward_odds_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    let m = prof.m;
    let l = dsq.len().saturating_sub(2);
    let om = OddsProfile::build(prof);

    // Two rolling rows of core-state odds (index 1..=m; 0 unused / stays 0).
    let mut pm = vec![0.0f32; m + 1];
    let mut pi = vec![0.0f32; m + 1];
    let mut pd = vec![0.0f32; m + 1];
    let mut cm = vec![0.0f32; m + 1];
    let mut ci = vec![0.0f32; m + 1];
    let mut cd = vec![0.0f32; m + 1];

    // Specials (odds). Row 0: N=1, B = N->B move, others 0.
    let (mut xn, mut xj, mut xc) = (1.0f32, 0.0f32, 0.0f32);
    let mut xb = om.n_move;
    let mut totscale = 0.0f64;

    for i in 1..=l {
        let me = &om.me[dsq[i] as usize];
        let mut xe = 0.0f32;

        // D is computed left-to-right this row (serial DD chain); M/I read the previous row.
        cd[1] = 0.0;
        for k in 1..=m {
            // Match: B->M and M/I/D(i-1, k-1) -> M(i, k), then * emission odds.
            let mm = xb * om.tbm[k - 1]
                + pm[k - 1] * om.tmm[k - 1]
                + pi[k - 1] * om.tim[k - 1]
                + pd[k - 1] * om.tdm[k - 1];
            cm[k] = mm * me[k];

            // Insert (emission odds 1.0): M/I(i-1, k) -> I(i, k).
            ci[k] = pm[k] * om.tmi[k] + pi[k] * om.tii[k];

            // Delete (same row): M/D(i, k-1) -> D(i, k).
            if k >= 2 {
                cd[k] = cm[k - 1] * om.tmd[k - 1] + cd[k - 1] * om.tdd[k - 1];
            }

            // E accumulates every Mk and Dk exit (esc = 1 odds in local mode).
            xe += cm[k] + cd[k];
        }

        // Specials (probability-space length model).
        xn *= om.n_loop;
        xc = xc * om.c_loop + xe * om.e_move;
        xj = xj * om.j_loop + xe * om.e_loop;
        xb = xj * om.j_move + xn * om.n_move;

        // Sparse rescaling: keep the row's dynamic range in check.
        if xe > RESCALE_THRESHOLD {
            let inv = 1.0 / xe;
            xn *= inv;
            xc *= inv;
            xj *= inv;
            xb *= inv;
            for k in 1..=m {
                cm[k] *= inv;
                ci[k] *= inv;
                cd[k] *= inv;
            }
            totscale += (xe as f64).ln();
        }

        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }

    // C->T close; flip the rescaled odds product back to nats.
    (totscale + (xc as f64 * om.c_move as f64).ln()) as f32
}

// ---------------------------------------------------------------------------------------------
// Stage 2: striped SSE kernel (`__m128`, 4 lanes) — a faithful port of HMMER's
// `forward_engine` (impl_sse/fwdback.c). Striping interleaves the M positions across lanes so
// the k-1 dependency becomes a cheap one-lane shift, letting the M/I/D recursion run 4-wide.
// ---------------------------------------------------------------------------------------------

/// Lane count of the SSE float vector.
#[cfg(target_arch = "x86_64")]
const VF: usize = 4;

/// Striped odds profile (`__m128` blocks), the layout HMMER's optimized profile (`om`) uses.
/// Position `p` (1..=M) maps to vector `q = (p-1) % Q`, lane `z = (p-1) / Q`.
///
/// The length-independent emission/transition blocks are built once via [`StripedProfile::new`];
/// only the N/C/J length-model odds change with target length, via [`StripedProfile::set_length`]
/// (cheap, no `exp`). This is what lets the windowing pipeline reuse one striped profile across
/// every tile instead of rebuilding it per window.
#[cfg(target_arch = "x86_64")]
#[derive(Clone)]
pub struct StripedProfile {
    q: usize,
    /// Match-emission odds per residue: `rfv[x]` is `Q` vectors.
    rfv: Vec<Vec<[f32; VF]>>,
    /// Interleaved transition odds: `Q` groups of 7 vectors, order [BM, MM, IM, DM, MD, MI, II].
    tfv: Vec<[f32; VF]>,
    /// D->D odds: `Q` vectors.
    dd: Vec<[f32; VF]>,
    n_loop: f32,
    n_move: f32,
    c_loop: f32,
    c_move: f32,
    j_loop: f32,
    j_move: f32,
    e_move: f32,
    e_loop: f32,
}

#[cfg(target_arch = "x86_64")]
impl StripedProfile {
    /// Build the striped odds profile (length-model specials initialized from `prof`'s current
    /// length config). Reuse it across windows; call [`StripedProfile::set_length`] per length.
    pub fn new(prof: &P7Profile) -> StripedProfile {
        Self::build(prof)
    }

    /// Reset the N/C/J length-model odds for target length `l` (HMMER's `p7_ReconfigLength`,
    /// in odds space — `ploop`/`pmove` are already probabilities, so no `exp` is needed).
    pub fn set_length(&mut self, l: usize, nj: f32) {
        let pmove = (2.0 + nj) / (l as f32 + 2.0 + nj);
        let ploop = 1.0 - pmove;
        self.n_loop = ploop;
        self.n_move = pmove;
        self.c_loop = ploop;
        self.c_move = pmove;
        self.j_loop = ploop;
        self.j_move = pmove;
    }

    /// Striped SSE Forward score in nats for `dsq` (1-based, sentinels at 0 and L+1).
    pub fn score_nats(&self, dsq: &[Dsq]) -> f32 {
        // SAFETY: SSE2 is part of the x86_64 baseline ABI, always present.
        unsafe { forward_striped_sse(self, dsq) }
    }

    fn build(prof: &P7Profile) -> StripedProfile {
        let m = prof.m;
        let q = m.div_ceil(VF);
        // Striped gather: lane z of vector qq holds node index `kb + z*Q`, odds = exp(tsc),
        // 0.0 when the index is out of range (HMMER's -inf padding).
        let gather = |ty: usize, kb_base: usize, valid_lt: usize| -> [f32; VF] {
            let mut v = [0.0f32; VF];
            for (z, slot) in v.iter_mut().enumerate() {
                let idx = kb_base + z * q;
                if idx < valid_lt {
                    *slot = prof.tsc(ty, idx).exp();
                }
            }
            v
        };

        let mut tfv = Vec::with_capacity(7 * q);
        for qq in 0..q {
            // BM/MM/IM/DM are rotated by -1 (kb = k-1 = qq); valid while kb+z*Q < M.
            tfv.push(gather(p7p::BM, qq, m));
            tfv.push(gather(p7p::MM, qq, m));
            tfv.push(gather(p7p::IM, qq, m));
            tfv.push(gather(p7p::DM, qq, m));
            // MD/MI/II use kb = k = qq+1; valid while kb+z*Q < M.
            tfv.push(gather(p7p::MD, qq + 1, m));
            tfv.push(gather(p7p::MI, qq + 1, m));
            tfv.push(gather(p7p::II, qq + 1, m));
        }
        let dd = (0..q).map(|qq| gather(p7p::DD, qq + 1, m)).collect();

        // Match emission: lane z holds position p = qq+1+z*Q (valid while p <= M).
        let rfv = prof
            .msc
            .iter()
            .map(|row| {
                (0..q)
                    .map(|qq| {
                        let mut v = [0.0f32; VF];
                        for (z, slot) in v.iter_mut().enumerate() {
                            let p = qq + 1 + z * q;
                            if p <= m {
                                *slot = row[p].exp();
                            }
                        }
                        v
                    })
                    .collect()
            })
            .collect();

        let xf = |st: usize, dir: usize| prof.xsc[st][dir].exp();
        StripedProfile {
            q,
            rfv,
            tfv,
            dd,
            n_loop: xf(xstate::N, xstate::LOOP),
            n_move: xf(xstate::N, xstate::MOVE),
            c_loop: xf(xstate::C, xstate::LOOP),
            c_move: xf(xstate::C, xstate::MOVE),
            j_loop: xf(xstate::J, xstate::LOOP),
            j_move: xf(xstate::J, xstate::MOVE),
            e_move: xf(xstate::E, xstate::MOVE),
            e_loop: xf(xstate::E, xstate::LOOP),
        }
    }
}

/// Striped SSE Forward score in nats. Mirrors [`forward_odds_nats`] (and thus
/// [`crate::fwd::forward_nats`]) to filter precision; the striped kernel is the SIMD path.
#[cfg(target_arch = "x86_64")]
pub fn forward_striped_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    StripedProfile::new(prof).score_nats(dsq)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[allow(clippy::needless_range_loop)] // striped DP indices (q, i) address several parallel rows
unsafe fn forward_striped_sse(sp: &StripedProfile, dsq: &[Dsq]) -> f32 {
    use std::arch::x86_64::*;

    let q = sp.q;
    let l = dsq.len().saturating_sub(2);
    let ld = |v: &[f32; VF]| _mm_loadu_ps(v.as_ptr());
    // Right-shift float vector by one lane, zero in: [a,b,c,d] -> [0,a,b,c] (esl_sse_rightshiftz).
    let rsz = |v: __m128| _mm_castsi128_ps(_mm_slli_si128::<4>(_mm_castps_si128(v)));

    let zero = _mm_setzero_ps();
    let mut pm = vec![zero; q];
    let mut pi = vec![zero; q];
    let mut pd = vec![zero; q];
    let mut cm = vec![zero; q];
    let mut ci = vec![zero; q];
    let mut cd = vec![zero; q];

    let (mut xn, mut xj, mut xc) = (1.0f32, 0.0f32, 0.0f32);
    let mut xb = sp.n_move;
    let mut totscale = 0.0f64;

    for i in 1..=l {
        let rp = &sp.rfv[dsq[i] as usize];
        let mut dcv = zero;
        let mut xev = zero;
        let xbv = _mm_set1_ps(xb);

        let mut mpv = rsz(pm[q - 1]);
        let mut dpv = rsz(pd[q - 1]);
        let mut ipv = rsz(pi[q - 1]);

        for qq in 0..q {
            let t = &sp.tfv[7 * qq..];
            // M(i,q) = (B->M + M->M + I->M + D->M) * emission.
            let mut sv = _mm_mul_ps(xbv, ld(&t[0]));
            sv = _mm_add_ps(sv, _mm_mul_ps(mpv, ld(&t[1])));
            sv = _mm_add_ps(sv, _mm_mul_ps(ipv, ld(&t[2])));
            sv = _mm_add_ps(sv, _mm_mul_ps(dpv, ld(&t[3])));
            sv = _mm_mul_ps(sv, ld(&rp[qq]));
            xev = _mm_add_ps(xev, sv);

            // Reload prev-row stripe qq (unshifted) for the next iteration's k-1.
            mpv = pm[qq];
            dpv = pd[qq];
            ipv = pi[qq];

            cm[qq] = sv;
            cd[qq] = dcv; // delayed D store (M->D only so far)
            dcv = _mm_mul_ps(sv, ld(&t[4])); // M->D into next stripe

            // I(i,q) = M(i-1,q)*M->I + I(i-1,q)*I->I (insert emission odds 1.0).
            ci[qq] = _mm_add_ps(_mm_mul_ps(mpv, ld(&t[5])), _mm_mul_ps(ipv, ld(&t[6])));
        }

        // DD propagation: one mandatory pass + full serialization (M < 100 ⇒ VF passes total).
        dcv = rsz(dcv);
        cd[0] = zero;
        for qq in 0..q {
            cd[qq] = _mm_add_ps(dcv, cd[qq]);
            dcv = _mm_mul_ps(cd[qq], ld(&sp.dd[qq]));
        }
        for _ in 1..VF {
            dcv = rsz(dcv);
            for qq in 0..q {
                cd[qq] = _mm_add_ps(dcv, cd[qq]);
                dcv = _mm_mul_ps(dcv, ld(&sp.dd[qq]));
            }
        }

        // E gets every Dk exit too.
        for qq in 0..q {
            xev = _mm_add_ps(xev, cd[qq]);
        }
        // Horizontal sum of xev -> xE.
        let h = _mm_add_ps(xev, _mm_shuffle_ps::<0x39>(xev, xev));
        let h = _mm_add_ps(h, _mm_shuffle_ps::<0x4e>(h, h));
        let xe = _mm_cvtss_f32(h);

        // Specials.
        xn *= sp.n_loop;
        xc = xc * sp.c_loop + xe * sp.e_move;
        xj = xj * sp.j_loop + xe * sp.e_loop;
        xb = xj * sp.j_move + xn * sp.n_move;

        // Sparse rescaling.
        if xe > RESCALE_THRESHOLD {
            let inv = 1.0 / xe;
            let iv = _mm_set1_ps(inv);
            xn *= inv;
            xc *= inv;
            xj *= inv;
            xb *= inv;
            for qq in 0..q {
                cm[qq] = _mm_mul_ps(cm[qq], iv);
                ci[qq] = _mm_mul_ps(ci[qq], iv);
                cd[qq] = _mm_mul_ps(cd[qq], iv);
            }
            totscale += (xe as f64).ln();
        }

        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }

    (totscale + (xc as f64 * sp.c_move as f64).ln()) as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fwd::forward_nats;
    use crate::profile::{bg_freqs, P7Profile};
    use easel_rs::Alphabet;

    /// The odds-space Forward must agree with the log-space DP to filter precision on a real
    /// model. Fixture: the embedded tRNA filter HMM + its consensus sequence.
    #[test]
    fn odds_forward_matches_log_forward() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../infernal-rs/tests/data/trna_fp7.hmm"
        );
        let hmm = crate::parse_p7_hmm(&std::fs::read_to_string(path).expect("read hmm"))
            .expect("parse hmm");
        let abc = Alphabet::rna();
        let bg = bg_freqs(hmm.k);

        // A few deterministic RNA sequences of differing length.
        let mut s: u64 = 0x51A7;
        let bases = [b'A', b'C', b'G', b'U'];
        for &len in &[20usize, 75, 160, 300] {
            let seq: String = (0..len)
                .map(|_| {
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                    bases[((s >> 40) & 3) as usize] as char
                })
                .collect();
            let dsq = abc.digitize(&seq).expect("digitize");
            let prof = P7Profile::config_local(&hmm, &abc, &bg, len);
            let log = forward_nats(&prof, &dsq);
            let odds = forward_odds_nats(&prof, &dsq);
            assert!(
                (log - odds).abs() < 0.02,
                "len {len}: log-space {log} vs odds-space {odds}"
            );
            #[cfg(target_arch = "x86_64")]
            {
                let striped = forward_striped_nats(&prof, &dsq);
                assert!(
                    (log - striped).abs() < 0.02,
                    "len {len}: log-space {log} vs striped-SSE {striped}"
                );
            }
        }
    }
}
