//! Striped-SIMD Viterbi filter — the max-plus counterpart of [`crate::fwd_simd`]'s striped
//! Forward. HMMER's real `p7_ViterbiFilter` (`impl_sse/vitfilter.c`) uses a reduced-precision
//! int16 kernel; here we keep the codebase's log-space **f32** representation (same as the
//! scalar [`crate::fwd::viterbi_nats`]) and stripe it, trading a little raw speed for exactness
//! and no quantization tuning. Because Viterbi is max-plus in log space it never underflows, so
//! — unlike the odds-space Forward — there is **no sparse rescaling**; the DP is pure add/max.
//!
//! The striping mirrors [`crate::fwd_simd::StripedProfile`] exactly (position `p` → vector
//! `q = (p-1) % Q`, lane `z = (p-1) / Q`), with two differences dictated by the semiring:
//!   * profile entries are stored as **log-odds** (`tsc`/`msc` directly, not `exp`), and
//!     out-of-range lanes are padded with `-inf` (the max-plus identity), not `0`;
//!   * the one-lane cross-stripe shift fills lane 0 with `-inf` instead of `0`.
//!
//! Verified against the scalar `viterbi_nats` oracle to f32 rounding (see tests).

use ornament_alphabet::Dsq;

use crate::profile::{p7p, xstate, P7Profile};

const NINF: f32 = f32::NEG_INFINITY;

/// Max SIMD lane width (AVX2 float). Each striped block is `[f32; VMAX]`; the active width
/// `vf` (4 = SSE, 8 = AVX2) selects how many lanes are populated and read.
#[cfg(target_arch = "x86_64")]
const VMAX: usize = 8;

/// Striped **log-odds** profile for the Viterbi filter (the max-plus analogue of
/// [`crate::fwd_simd::StripedProfile`]). Length-independent emission/transition blocks are
/// built once via [`StripedVitProfile::new`]; only the N/C/J length-model log-odds change per
/// target length via [`StripedVitProfile::set_length`], so the windowing pipeline reuses one
/// profile across every tile.
#[cfg(target_arch = "x86_64")]
#[derive(Clone)]
pub struct StripedVitProfile {
    /// Active vector width: 8 (AVX2) or 4 (SSE). Only lanes `0..vf` of each block are used.
    vf: usize,
    q: usize,
    /// Match-emission log-odds per residue: `rfv[x]` is `Q` vectors (`-inf` padded).
    rfv: Vec<Vec<[f32; VMAX]>>,
    /// Interleaved transition log-odds: `Q` groups of 7 vectors [BM, MM, IM, DM, MD, MI, II].
    tfv: Vec<[f32; VMAX]>,
    /// D->D log-odds: `Q` vectors.
    dd: Vec<[f32; VMAX]>,
    // Special-state log-odds (length model in log space).
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
impl StripedVitProfile {
    /// Build the striped log-odds profile (length-model specials from `prof`'s current length
    /// config), striped for the widest SIMD the CPU supports (AVX2 → 8 lanes, else SSE → 4).
    pub fn new(prof: &P7Profile) -> StripedVitProfile {
        let vf = if std::is_x86_feature_detected!("avx2") {
            8
        } else {
            4
        };
        Self::build(prof, vf)
    }

    /// Reset the N/C/J length-model log-odds for target length `l` (HMMER's `p7_ReconfigLength`),
    /// mirroring [`crate::fwd_simd::StripedProfile::set_length`] but keeping the values in log
    /// space (`ln` of the loop/move probabilities).
    pub fn set_length(&mut self, l: usize, nj: f32) {
        let pmove = (2.0 + nj) / (l as f32 + 2.0 + nj);
        let ploop = 1.0 - pmove;
        let (lln, lmv) = (ploop.ln(), pmove.ln());
        self.n_loop = lln;
        self.n_move = lmv;
        self.c_loop = lln;
        self.c_move = lmv;
        self.j_loop = lln;
        self.j_move = lmv;
    }

    /// Striped Viterbi score in nats for `dsq` (1-based, sentinels at 0 and L+1). Dispatches to
    /// the AVX2 (8-lane) or SSE (4-lane) kernel matching how the profile was striped. Matches the
    /// scalar [`crate::fwd::viterbi_nats`] to f32 rounding.
    pub fn score_nats(&self, dsq: &[Dsq]) -> f32 {
        // SAFETY: the width was chosen from `is_x86_feature_detected!` at build time; SSE2 is
        // baseline and AVX2 was confirmed present when `vf == 8`.
        unsafe {
            if self.vf == 8 {
                viterbi_striped_avx2(self, dsq)
            } else {
                viterbi_striped_sse(self, dsq)
            }
        }
    }

    fn build(prof: &P7Profile, vf: usize) -> StripedVitProfile {
        let m = prof.m;
        let q = m.div_ceil(vf);
        // Striped gather: lane z (0..vf) of vector qq holds node index `kb_base + z*Q`, value =
        // log-odds `tsc`, `-inf` when out of range (max-plus identity). Lanes vf..VMAX stay -inf.
        let gather = |ty: usize, kb_base: usize, valid_lt: usize| -> [f32; VMAX] {
            let mut v = [NINF; VMAX];
            for (z, slot) in v.iter_mut().enumerate().take(vf) {
                let idx = kb_base + z * q;
                if idx < valid_lt {
                    *slot = prof.tsc(ty, idx);
                }
            }
            v
        };

        let mut tfv = Vec::with_capacity(7 * q);
        for qq in 0..q {
            // BM/MM/IM/DM are indexed by k-1 = qq; valid while kb+z*Q < M.
            tfv.push(gather(p7p::BM, qq, m));
            tfv.push(gather(p7p::MM, qq, m));
            tfv.push(gather(p7p::IM, qq, m));
            tfv.push(gather(p7p::DM, qq, m));
            // MD/MI/II are indexed by k = qq+1; valid while kb+z*Q < M.
            tfv.push(gather(p7p::MD, qq + 1, m));
            tfv.push(gather(p7p::MI, qq + 1, m));
            tfv.push(gather(p7p::II, qq + 1, m));
        }
        let dd = (0..q).map(|qq| gather(p7p::DD, qq + 1, m)).collect();

        // Match emission log-odds: lane z holds position p = qq+1+z*Q (valid while p <= M).
        let rfv = prof
            .msc
            .iter()
            .map(|row| {
                (0..q)
                    .map(|qq| {
                        let mut v = [NINF; VMAX];
                        for (z, slot) in v.iter_mut().enumerate().take(vf) {
                            let p = qq + 1 + z * q;
                            if p <= m {
                                *slot = row[p];
                            }
                        }
                        v
                    })
                    .collect()
            })
            .collect();

        let xf = |st: usize, dir: usize| prof.xsc[st][dir];
        StripedVitProfile {
            vf,
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

/// Convenience: striped Viterbi score in nats (builds a fresh profile). Mirrors
/// [`crate::fwd_simd::forward_striped_nats`].
#[cfg(target_arch = "x86_64")]
pub fn viterbi_striped_nats(prof: &P7Profile, dsq: &[Dsq]) -> f32 {
    StripedVitProfile::new(prof).score_nats(dsq)
}

/// Striped SSE Viterbi (`__m128`, 4 lanes). Max-plus port of [`crate::fwd_simd::forward_striped_sse`]:
/// `*`→`+` (transition/emission compose additively) and `+`→`max` (alternative sub-parses).
/// The cross-stripe shift fills lane 0 with `-inf`, and there is no sparse rescaling.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[allow(clippy::needless_range_loop)]
unsafe fn viterbi_striped_sse(sp: &StripedVitProfile, dsq: &[Dsq]) -> f32 {
    use std::arch::x86_64::*;

    let q = sp.q;
    let l = dsq.len().saturating_sub(2);
    let ld = |v: &[f32; VMAX]| _mm_loadu_ps(v.as_ptr());
    let ninf = _mm_set1_ps(NINF);
    // -inf only in lane 0, +0 elsewhere; OR-ing into a slli-zeroed lane 0 sets it to -inf
    // (0-bits OR NINF-bits = NINF-bits) while leaving the shifted lanes untouched. sse2-only.
    let ninf0 = _mm_set_ps(0.0, 0.0, 0.0, NINF);
    // Right-shift by one lane (zero-fill), then patch lane 0 to -inf: [a,b,c,d] -> [-inf,a,b,c].
    let rsz = |v: __m128| {
        _mm_or_ps(
            _mm_castsi128_ps(_mm_slli_si128::<4>(_mm_castps_si128(v))),
            ninf0,
        )
    };

    let mut pm = vec![ninf; q];
    let mut pi = vec![ninf; q];
    let mut pd = vec![ninf; q];
    let mut cm = vec![ninf; q];
    let mut ci = vec![ninf; q];
    let mut cd = vec![ninf; q];

    let (mut xn, mut xj, mut xc) = (0.0f32, NINF, NINF);
    let mut xb = sp.n_move;

    for i in 1..=l {
        let rp = &sp.rfv[dsq[i] as usize];
        let mut dcv = ninf;
        let mut xev = ninf;
        let xbv = _mm_set1_ps(xb);

        let mut mpv = rsz(pm[q - 1]);
        let mut dpv = rsz(pd[q - 1]);
        let mut ipv = rsz(pi[q - 1]);

        for qq in 0..q {
            let t = &sp.tfv[7 * qq..];
            // M(i,q) = max(B->M, M->M, I->M, D->M) + emission.
            let mut sv = _mm_add_ps(xbv, ld(&t[0]));
            sv = _mm_max_ps(sv, _mm_add_ps(mpv, ld(&t[1])));
            sv = _mm_max_ps(sv, _mm_add_ps(ipv, ld(&t[2])));
            sv = _mm_max_ps(sv, _mm_add_ps(dpv, ld(&t[3])));
            sv = _mm_add_ps(sv, ld(&rp[qq]));
            xev = _mm_max_ps(xev, sv);

            // Reload prev-row stripe qq (unshifted) for the next iteration's k-1.
            mpv = pm[qq];
            dpv = pd[qq];
            ipv = pi[qq];

            cm[qq] = sv;
            cd[qq] = dcv; // delayed D store (M->D only so far)
            dcv = _mm_add_ps(sv, ld(&t[4])); // M->D into next stripe

            // I(i,q) = max(M(i-1,q)+M->I, I(i-1,q)+I->I) (insert emission log-odds 0).
            ci[qq] = _mm_max_ps(_mm_add_ps(mpv, ld(&t[5])), _mm_add_ps(ipv, ld(&t[6])));
        }

        // DD propagation: max-plus, full serialization (4 passes for SSE, M < 100).
        dcv = rsz(dcv);
        cd[0] = ninf;
        for qq in 0..q {
            cd[qq] = _mm_max_ps(dcv, cd[qq]);
            dcv = _mm_add_ps(cd[qq], ld(&sp.dd[qq]));
        }
        for _ in 1..4 {
            dcv = rsz(dcv);
            for qq in 0..q {
                cd[qq] = _mm_max_ps(dcv, cd[qq]);
                dcv = _mm_add_ps(dcv, ld(&sp.dd[qq]));
            }
        }

        // E gets every Dk exit too.
        for qq in 0..q {
            xev = _mm_max_ps(xev, cd[qq]);
        }
        // Horizontal max of xev -> xE.
        let h = _mm_max_ps(xev, _mm_shuffle_ps::<0x39>(xev, xev));
        let h = _mm_max_ps(h, _mm_shuffle_ps::<0x4e>(h, h));
        let xe = _mm_cvtss_f32(h);

        // Specials (log-space max-plus length model).
        xn += sp.n_loop;
        xc = (xc + sp.c_loop).max(xe + sp.e_move);
        xj = (xj + sp.j_loop).max(xe + sp.e_loop);
        xb = (xn + sp.n_move).max(xj + sp.j_move);

        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }

    xc + sp.c_move
}

/// AVX2 (8-lane `__m256`) striped Viterbi — same algorithm as [`viterbi_striped_sse`] at double
/// the width. Requires the profile to be striped for `vf = 8`.
///
/// # Safety
/// Requires the `avx2` target feature (the build picked `vf == 8` only after detecting it).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::needless_range_loop)]
unsafe fn viterbi_striped_avx2(sp: &StripedVitProfile, dsq: &[Dsq]) -> f32 {
    use std::arch::x86_64::*;

    let q = sp.q;
    let l = dsq.len().saturating_sub(2);
    let ld = |v: &[f32; VMAX]| _mm256_loadu_ps(v.as_ptr());
    let ninf = _mm256_set1_ps(NINF);
    // Rotate by one lane then set lane 0 to -inf (blend, since AVX2 has no cross-128 byte shift).
    let shuf = _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6);
    let rsz = |v: __m256| _mm256_blend_ps::<0x01>(_mm256_permutevar8x32_ps(v, shuf), ninf);

    let mut pm = vec![ninf; q];
    let mut pi = vec![ninf; q];
    let mut pd = vec![ninf; q];
    let mut cm = vec![ninf; q];
    let mut ci = vec![ninf; q];
    let mut cd = vec![ninf; q];

    let (mut xn, mut xj, mut xc) = (0.0f32, NINF, NINF);
    let mut xb = sp.n_move;

    for i in 1..=l {
        let rp = &sp.rfv[dsq[i] as usize];
        let mut dcv = ninf;
        let mut xev = ninf;
        let xbv = _mm256_set1_ps(xb);

        let mut mpv = rsz(pm[q - 1]);
        let mut dpv = rsz(pd[q - 1]);
        let mut ipv = rsz(pi[q - 1]);

        for qq in 0..q {
            let t = &sp.tfv[7 * qq..];
            let mut sv = _mm256_add_ps(xbv, ld(&t[0]));
            sv = _mm256_max_ps(sv, _mm256_add_ps(mpv, ld(&t[1])));
            sv = _mm256_max_ps(sv, _mm256_add_ps(ipv, ld(&t[2])));
            sv = _mm256_max_ps(sv, _mm256_add_ps(dpv, ld(&t[3])));
            sv = _mm256_add_ps(sv, ld(&rp[qq]));
            xev = _mm256_max_ps(xev, sv);

            mpv = pm[qq];
            dpv = pd[qq];
            ipv = pi[qq];

            cm[qq] = sv;
            cd[qq] = dcv;
            dcv = _mm256_add_ps(sv, ld(&t[4]));

            ci[qq] = _mm256_max_ps(_mm256_add_ps(mpv, ld(&t[5])), _mm256_add_ps(ipv, ld(&t[6])));
        }

        // DD propagation: 8 passes total for AVX2 (full serialization, M < 100).
        dcv = rsz(dcv);
        cd[0] = ninf;
        for qq in 0..q {
            cd[qq] = _mm256_max_ps(dcv, cd[qq]);
            dcv = _mm256_add_ps(cd[qq], ld(&sp.dd[qq]));
        }
        for _ in 1..8 {
            dcv = rsz(dcv);
            for qq in 0..q {
                cd[qq] = _mm256_max_ps(dcv, cd[qq]);
                dcv = _mm256_add_ps(dcv, ld(&sp.dd[qq]));
            }
        }

        for qq in 0..q {
            xev = _mm256_max_ps(xev, cd[qq]);
        }
        // Horizontal max of 8 lanes -> xE.
        let lo = _mm256_castps256_ps128(xev);
        let hi = _mm256_extractf128_ps::<1>(xev);
        let s = _mm_max_ps(lo, hi);
        let s = _mm_max_ps(s, _mm_shuffle_ps::<0x39>(s, s));
        let s = _mm_max_ps(s, _mm_shuffle_ps::<0x4e>(s, s));
        let xe = _mm_cvtss_f32(s);

        xn += sp.n_loop;
        xc = (xc + sp.c_loop).max(xe + sp.e_move);
        xj = (xj + sp.j_loop).max(xe + sp.e_loop);
        xb = (xn + sp.n_move).max(xj + sp.j_move);

        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }

    xc + sp.c_move
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;
    use crate::fwd::viterbi_nats;
    use crate::profile::{bg_freqs, P7Profile};
    use ornament_alphabet::Alphabet;

    /// The striped Viterbi must agree with the scalar max-plus DP to f32 rounding on a real
    /// model. Fixture: the embedded tRNA filter HMM. (Max-plus has no rescaling, so agreement
    /// is far tighter than Forward's odds-space tolerance.)
    #[test]
    fn striped_viterbi_matches_scalar() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../ornament-scfg/tests/data/trna_fp7.hmm"
        );
        let hmm = crate::parse_p7_hmm(&std::fs::read_to_string(path).expect("read hmm"))
            .expect("parse hmm");
        let abc = Alphabet::rna();
        let bg = bg_freqs(hmm.k);

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
            let scalar = viterbi_nats(&prof, &dsq);

            let widths: &[usize] = if std::is_x86_feature_detected!("avx2") {
                &[4, 8]
            } else {
                &[4]
            };
            for &vf in widths {
                let striped = StripedVitProfile::build(&prof, vf).score_nats(&dsq);
                assert!(
                    (scalar - striped).abs() < 0.01,
                    "len {len} vf {vf}: scalar {scalar} vs striped {striped}"
                );
            }
        }
    }
}
