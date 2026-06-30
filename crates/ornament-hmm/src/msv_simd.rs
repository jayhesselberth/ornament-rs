//! Striped reduced-precision (uint8) MSV filter — a port of HMMER's `p7_MSVFilter`
//! (`impl_sse/msvfilter.c`), the cheapest-first filter in the cascade.
//!
//! MSV is ungapped multi-segment Viterbi (no I/D states, flat B→Mk entry). Because it is
//! max-plus, it fits unsigned-saturating 8-bit arithmetic directly — **16 lanes** per SSE
//! vector (4× the f32 Forward), with a much lighter recurrence. Scores live in offset
//! "third-bit" units (`scale_b = 3/ln2`, `base_b = 190`): `0..255` maps to roughly
//! `−63..+22` bits, so a strong hit *saturates*. Saturation is detected and reported as
//! `+∞` so the caller treats the tile as a pass and lets a higher-precision filter (Forward)
//! score it — exactly HMMER's design. Junk tiles score finite and are pruned cheaply.
//!
//! The recovered raw score (`(xJ − tjb_b − base_b)/scale_b − 3.0`, nats) is null-corrected to
//! bits as `(raw − null_one(L))/ln2` — the same convention as the scalar [`crate::msv_bits`], so
//! it's on the scale the calibrated [`crate::pvalue::msv_pvalue`] expects. (The uint8 floor makes
//! it diverge from the float MSV on *weak* tiles, but those are pruned anyway; near the filter
//! threshold the two agree to within quantization, which the cascade's `F1` margin absorbs.)

#[cfg(target_arch = "x86_64")]
use crate::profile::P7Profile;
#[cfg(target_arch = "x86_64")]
use ornament_alphabet::Dsq;

/// Max uint8 lane count (AVX2). Each striped block is stored `[u8; VBMAX]`; the active width
/// `vb` (16 = SSE, 32 = AVX2) selects how many lanes are populated and read.
#[cfg(target_arch = "x86_64")]
const VBMAX: usize = 32;

/// `unbiased_byteify`: round a (log-prob) score to a non-negative uint8 cost. `255` = −∞.
#[cfg(target_arch = "x86_64")]
fn unbiased_byteify(scale_b: f32, sc: f32) -> u8 {
    let cost = -(scale_b * sc).round();
    if cost > 255.0 {
        255
    } else {
        cost.max(0.0) as u8
    }
}

/// `biased_byteify`: emission score → biased uint8 cost (`bias_b − round(scale_b·sc)`), with
/// HMMER's two's-complement wraparound (`(uint8)cost + bias_b`) and `255` = −∞. In the DP we
/// add `bias_b` then subtract this, recovering `round(scale_b·sc)` in saturating arithmetic.
#[cfg(target_arch = "x86_64")]
fn biased_byteify(scale_b: f32, bias_b: u8, sc: f32) -> u8 {
    let cost = -(scale_b * sc).round();
    if cost > 255.0 - bias_b as f32 {
        255
    } else {
        (cost as i32 as u8).wrapping_add(bias_b)
    }
}

/// Striped uint8 MSV profile. Built once from a configured [`P7Profile`] (reusing its match
/// log-odds), reused across windows; [`MsvProfile::set_length`] resets only the length-model
/// byte `tjb_b`.
#[cfg(target_arch = "x86_64")]
#[derive(Clone)]
pub struct MsvProfile {
    /// Active lane width: 32 (AVX2) or 16 (SSE). Only lanes `0..vb` of each block are used.
    vb: usize,
    q: usize,
    scale_b: f32,
    base_b: u8,
    bias_b: u8,
    tbm_b: u8, // B->Mk entry (from M)
    tec_b: u8, // E->C / E->J (constant, multihit)
    tjb_b: u8, // N/C/J move (length-dependent)
    /// `rbv[x]` = `Q` striped blocks of biased emission costs for residue `x` (lanes `vb..` = 255).
    rbv: Vec<Vec<[u8; VBMAX]>>,
}

#[cfg(target_arch = "x86_64")]
impl MsvProfile {
    /// Build the uint8 MSV profile, striped for the widest SIMD that actually helps. AVX2's
    /// 32-lane kernel only wins once the model is long enough to amortize its pricier cross-lane
    /// shift — at small `M` (few stripes, e.g. tRNA's M≈71) it's marginally *slower* than SSE
    /// (measured), so we keep SSE 16-lane below `M = 128` and switch to AVX2 above it. Length
    /// model from `prof.l`.
    pub fn new(prof: &P7Profile) -> MsvProfile {
        let vb = if prof.m >= 128 && std::is_x86_feature_detected!("avx2") {
            32
        } else {
            16
        };
        Self::build(prof, vb)
    }

    /// Build forcing a specific SIMD lane width (16 = SSE, 32 = AVX2). For benchmarks/tests;
    /// production code uses [`MsvProfile::new`] (auto-selects the widest supported). Scoring a
    /// `vb = 32` profile on a non-AVX2 CPU is unsound, so only use 32 when AVX2 is present.
    pub fn with_lanes(prof: &P7Profile, lanes: usize) -> MsvProfile {
        Self::build(prof, lanes)
    }

    fn build(prof: &P7Profile, vb: usize) -> MsvProfile {
        let m = prof.m;
        let kp = prof.kp;
        let q = m.div_ceil(vb);
        let scale_b = 3.0 / std::f32::consts::LN_2;
        let base_b = 190u8;

        // bias_b from the maximum match-emission log-odds (over finite cells).
        let mut max = 0.0f32;
        for row in prof.msc.iter() {
            for &s in row[1..=m].iter() {
                if s.is_finite() {
                    max = max.max(s);
                }
            }
        }
        let bias_b = unbiased_byteify(scale_b, -max);

        // Striped biased emission bytes: lane z (0..vb) of vector q holds position p = q+1+z*Q.
        // Out-of-range positions and the unused `vb..VBMAX` lanes stay 255 (the -inf cost).
        let rbv: Vec<Vec<[u8; VBMAX]>> = (0..kp)
            .map(|x| {
                (0..q)
                    .map(|qq| {
                        let mut v = [255u8; VBMAX];
                        for (z, b) in v.iter_mut().enumerate().take(vb) {
                            let p = qq + 1 + z * q;
                            if p <= m {
                                *b = biased_byteify(scale_b, bias_b, prof.msc[x][p]);
                            }
                        }
                        v
                    })
                    .collect()
            })
            .collect();

        let tbm_b = unbiased_byteify(scale_b, (2.0 / (m as f32 * (m as f32 + 1.0))).ln());
        let tec_b = unbiased_byteify(scale_b, 0.5f32.ln());
        let mut sp = MsvProfile {
            vb,
            q,
            scale_b,
            base_b,
            bias_b,
            tbm_b,
            tec_b,
            tjb_b: 0,
            rbv,
        };
        sp.set_length(prof.l);
        sp
    }

    /// Reset the N/C/J length-model byte for target length `l` (`tjb_b = byteify(log 3/(L+3))`).
    pub fn set_length(&mut self, l: usize) {
        self.tjb_b = unbiased_byteify(self.scale_b, (3.0 / (l as f32 + 3.0)).ln());
    }

    /// MSV **bit** score for `dsq` (1-based, sentinels at 0 and L+1). Dispatches to the AVX2
    /// (32-lane) or SSE (16-lane) kernel matching how the profile was striped. Returns `+∞` on
    /// uint8 overflow (a saturating-strong hit) — the caller treats that as a pass.
    pub fn score_bits(&self, dsq: &[Dsq]) -> f32 {
        // SAFETY: width was chosen from `is_x86_feature_detected!` at build time; SSE2 is
        // baseline and AVX2 was confirmed present when `vb == 32`.
        unsafe {
            if self.vb == 32 {
                msv_avx2(self, dsq)
            } else {
                msv_sse(self, dsq)
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[allow(clippy::needless_range_loop)]
unsafe fn msv_sse(sp: &MsvProfile, dsq: &[Dsq]) -> f32 {
    use std::arch::x86_64::*;

    let q = sp.q;
    let l = dsq.len().saturating_sub(2);
    let biasv = _mm_set1_epi8(sp.bias_b as i8);
    let basev = _mm_set1_epi8(sp.base_b as i8);
    let tjbmv = _mm_set1_epi8(sp.tjb_b.wrapping_add(sp.tbm_b) as i8);
    let tecv = _mm_set1_epi8(sp.tec_b as i8);

    let mut dp = vec![_mm_setzero_si128(); q];
    let mut xjv = _mm_setzero_si128(); // xJ = 0
    let mut xbv = _mm_subs_epu8(basev, tjbmv); // xB = base - (tjb+tbm)

    for i in 1..=l {
        let rsc = &sp.rbv[dsq[i] as usize];
        let mut xev = _mm_setzero_si128();
        // Right-shift by one byte (striped k-1 wrap; zeros = -inf shift on).
        let mut mpv = _mm_slli_si128::<1>(dp[q - 1]);

        for qq in 0..q {
            // sv = max(M(i-1,k-1), xB) + emission   (offset arithmetic: +bias then -cost)
            let mut sv = _mm_max_epu8(mpv, xbv);
            sv = _mm_adds_epu8(sv, biasv);
            sv = _mm_subs_epu8(sv, _mm_loadu_si128(rsc[qq].as_ptr() as *const __m128i));
            xev = _mm_max_epu8(xev, sv);
            mpv = dp[qq];
            dp[qq] = sv;
        }

        // Horizontal max over the 16 bytes (scalar reduce — cheap, once per row).
        let mut bytes = [0u8; 16];
        _mm_storeu_si128(bytes.as_mut_ptr() as *mut __m128i, xev);
        let xe = bytes.iter().copied().max().unwrap();

        // Overflow: a lane saturated (xE + bias == 255) ⇒ strong hit ⇒ pass it on (+inf).
        if xe.saturating_add(sp.bias_b) == 255 {
            return f32::INFINITY;
        }

        // Specials: xJ = max(xJ, xE - tec); xB = max(base, xJ) - (tjb+tbm).
        let xev_b = _mm_set1_epi8(xe as i8);
        xjv = _mm_max_epu8(xjv, _mm_subs_epu8(xev_b, tecv));
        xbv = _mm_subs_epu8(_mm_max_epu8(basev, xjv), tjbmv);
    }

    // Recover the raw MSV score in nats (HMMER's MSVFilter recovery), then null-correct to bits
    // exactly as the scalar `msv_bits` does, so the result is on the scale `msv_pvalue` expects.
    let xj = (_mm_extract_epi16::<0>(xjv) & 0xff) as u8;
    recover_bits(sp, xj, l)
}

/// Shared score recovery (SSE/AVX2): byte → null-corrected bits on the `msv_pvalue` scale.
#[cfg(target_arch = "x86_64")]
#[inline]
fn recover_bits(sp: &MsvProfile, xj: u8, l: usize) -> f32 {
    let raw_nats = ((xj as f32 - sp.tjb_b as f32) - sp.base_b as f32) / sp.scale_b - 3.0;
    (raw_nats - crate::profile::null_one(l)) / std::f32::consts::LN_2
}

/// AVX2 (32-lane) uint8 MSV — the same algorithm as [`msv_sse`] at double the width. Requires
/// the profile striped for `vb = 32`.
///
/// # Safety
/// Requires the `avx2` target feature (the build picked `vb == 32` only after detecting it).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::needless_range_loop)]
unsafe fn msv_avx2(sp: &MsvProfile, dsq: &[Dsq]) -> f32 {
    use std::arch::x86_64::*;

    let q = sp.q;
    let l = dsq.len().saturating_sub(2);
    let biasv = _mm256_set1_epi8(sp.bias_b as i8);
    let basev = _mm256_set1_epi8(sp.base_b as i8);
    let tjbmv = _mm256_set1_epi8(sp.tjb_b.wrapping_add(sp.tbm_b) as i8);
    let tecv = _mm256_set1_epi8(sp.tec_b as i8);

    // Shift the whole 256-bit register left by one byte, zero in (the striped k-1 wrap). AVX2's
    // `slli_si256` only shifts within each 128-bit lane, so bring lane0's top byte across into
    // lane1's low byte via `permute2x128` (→ [lane0 | 0]) + a per-lane `alignr`.
    let shift1 = |v: __m256i| -> __m256i {
        let t = _mm256_permute2x128_si256::<0x08>(v, v); // [hi = lane0(v), lo = 0]
        _mm256_alignr_epi8::<15>(v, t)
    };

    let mut dp = vec![_mm256_setzero_si256(); q];
    let mut xjv = _mm256_setzero_si256();
    let mut xbv = _mm256_subs_epu8(basev, tjbmv);

    for i in 1..=l {
        let rsc = &sp.rbv[dsq[i] as usize];
        let mut xev = _mm256_setzero_si256();
        let mut mpv = shift1(dp[q - 1]);

        for qq in 0..q {
            let mut sv = _mm256_max_epu8(mpv, xbv);
            sv = _mm256_adds_epu8(sv, biasv);
            sv = _mm256_subs_epu8(sv, _mm256_loadu_si256(rsc[qq].as_ptr() as *const __m256i));
            xev = _mm256_max_epu8(xev, sv);
            mpv = dp[qq];
            dp[qq] = sv;
        }

        let mut bytes = [0u8; 32];
        _mm256_storeu_si256(bytes.as_mut_ptr() as *mut __m256i, xev);
        let xe = bytes.iter().copied().max().unwrap();

        if xe.saturating_add(sp.bias_b) == 255 {
            return f32::INFINITY;
        }

        let xev_b = _mm256_set1_epi8(xe as i8);
        xjv = _mm256_max_epu8(xjv, _mm256_subs_epu8(xev_b, tecv));
        xbv = _mm256_subs_epu8(_mm256_max_epu8(basev, xjv), tjbmv);
    }

    let xj = (_mm256_extract_epi8::<0>(xjv) & 0xff) as u8;
    recover_bits(sp, xj, l)
}

#[cfg(test)]
#[cfg(target_arch = "x86_64")]
mod tests {
    use super::*;
    use crate::msv_bits;
    use crate::profile::{bg_freqs, P7Profile};
    use ornament_alphabet::Alphabet;

    /// Behavioural validation (uint8 MSV is HMMER's reduced-precision filter — it floors weak
    /// scores, so it does NOT equal the float MSV bit-for-bit; it's the filter `msv_pvalue` is
    /// calibrated to). What must hold: a strong consensus hit passes the MSV threshold (so the
    /// cascade never drops it), and junk scores far lower (so it prunes). End-to-end safety is
    /// covered by the pipeline parity test.
    #[test]
    fn int8_msv_passes_hits_and_separates_junk() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../ornament-scfg/tests/data/trna_fp7.hmm"
        );
        let hmm = crate::parse_p7_hmm(&std::fs::read_to_string(path).expect("read hmm"))
            .expect("parse hmm");
        let abc = Alphabet::rna();
        let bg = bg_freqs(hmm.k);

        // Strong consensus tRNA: must pass MSV F1 (overflow → +inf, or a high finite score) so
        // the cascade always hands it to Forward.
        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../ornament-scfg/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        let dcons = abc.digitize(&cons).expect("digitize");
        let pcons = P7Profile::config_local(&hmm, &abc, &bg, cons.len());
        let cons_msv = MsvProfile::new(&pcons).score_bits(&dcons);
        assert!(
            crate::msv_pvalue(&hmm, cons_msv) <= 0.02,
            "consensus must pass MSV F1 (got bits {cons_msv}, scalar {})",
            msv_bits(&hmm, &abc, &dcons)
        );

        // A signal-free random tile must score far below the hit and fail F1 (gets pruned).
        let mut s: u64 = 0xA5A5;
        let bases = [b'A', b'C', b'G', b'U'];
        let seq: String = (0..436)
            .map(|_| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                bases[((s >> 40) & 3) as usize] as char
            })
            .collect();
        let drand = abc.digitize(&seq).expect("digitize");
        let prand = P7Profile::config_local(&hmm, &abc, &bg, 436);
        let rand_msv = MsvProfile::new(&prand).score_bits(&drand);
        assert!(
            rand_msv.is_finite(),
            "junk tile shouldn't spuriously overflow"
        );
        assert!(
            crate::msv_pvalue(&hmm, rand_msv) > 0.02,
            "junk tile should be pruned by MSV F1 (got bits {rand_msv})"
        );
        assert!(
            cons_msv > rand_msv,
            "hit {cons_msv} must outscore junk {rand_msv}"
        );
    }

    /// SSE (16-lane) and AVX2 (32-lane) must produce the *same* score — the striping width only
    /// changes the traversal order of an associative max-plus DP, not the result.
    #[test]
    fn msv_sse_and_avx2_agree() {
        let hmm = crate::parse_p7_hmm(
            &std::fs::read_to_string(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../ornament-scfg/tests/data/trna_fp7.hmm"
            ))
            .expect("read hmm"),
        )
        .expect("parse hmm");
        let abc = Alphabet::rna();
        let bg = bg_freqs(hmm.k);
        if !std::is_x86_feature_detected!("avx2") {
            return; // no AVX2 to compare against
        }
        let mut s: u64 = 0xBEEF;
        let bases = [b'A', b'C', b'G', b'U'];
        for &len in &[20usize, 73, 200, 436, 900] {
            let seq: String = (0..len)
                .map(|_| {
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                    bases[((s >> 40) & 3) as usize] as char
                })
                .collect();
            let dsq = abc.digitize(&seq).expect("digitize");
            let prof = P7Profile::config_local(&hmm, &abc, &bg, len);
            let sse = MsvProfile::build(&prof, 16).score_bits(&dsq);
            let avx2 = MsvProfile::build(&prof, 32).score_bits(&dsq);
            assert_eq!(
                sse.to_bits(),
                avx2.to_bits(),
                "len {len}: SSE {sse} vs AVX2 {avx2}"
            );
        }
    }
}
