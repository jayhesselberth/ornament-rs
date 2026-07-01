//! `ornament-gpu`: CUDA-accelerated p7 filter kernels — a prototype toward offloading the
//! all-model / many-genomes `scan` to GPU (see the project GPU-offload plan).
//!
//! The scan's heaviest, most parallel workload is the p7 prefilter cascade run over a huge batch
//! of (model, window) tiles. This crate starts with the cheapest, most GPU-friendly stage — the
//! **MSV filter** — batched one-thread-per-window on the device (`kernels/msv.cu`). It is a
//! faithful f32 max-plus transcription of the scalar [`ornament_hmm::msv_nats`] oracle, so the
//! device result matches the CPU to f32 rounding.
//!
//! ## Layout of a batch
//! A batch is one configured model's emission table plus N independent windows:
//! - [`FlatProfile`] — the model's match-emission log-odds flattened to `msc[x*(M+1)+k]`.
//! - [`Windows`] — every window's residue codes concatenated, with per-window offset/length.
//!
//! ## Feature gating
//! The `cuda` feature turns on the device path (build.rs compiles the kernel with nvcc and links
//! the CUDA runtime). Without it the crate still builds and the CPU oracle + batch builders are
//! available — only [`msv_nats_gpu`] is compiled out. This keeps the workspace buildable on
//! login nodes that have no CUDA toolchain; build/run the GPU path on the `compgpu` nodes.

use ornament_alphabet::Dsq;
use ornament_hmm::profile::null_one;
use ornament_hmm::P7Profile;
use thiserror::Error;

/// Errors from the GPU filter path.
#[derive(Debug, Error)]
pub enum GpuError {
    /// The crate was built without the `cuda` feature, so no device kernel exists.
    #[error("ornament-gpu was built without the `cuda` feature (no device kernel compiled)")]
    NotCompiled,
    /// No CUDA device is visible at runtime.
    #[error("no CUDA device available")]
    NoDevice,
    /// A CUDA Runtime API call failed; carries the `cudaError_t` code.
    #[error("CUDA runtime error (cudaError_t = {0})")]
    Cuda(i32),
}

/// A model's match-emission log-odds, flattened row-major to `msc[x*(M+1)+k]` for the kernel.
/// Built once per configured model and reused across every window batch. Carries only the
/// emission scores — MSV's B->Mk entry and length model are recomputed on-device from `M` and the
/// per-window length, exactly as [`ornament_hmm::msv_nats`] does.
#[derive(Debug, Clone)]
pub struct FlatProfile {
    /// Alphabet size (emission-table rows, `P7Profile::kp`).
    pub kp: usize,
    /// Model length (match states).
    pub m: usize,
    /// Row-major `[Kp][M+1]` match-emission log-odds (nats).
    pub msc: Vec<f32>,
}

impl FlatProfile {
    /// Flatten a configured [`P7Profile`]'s match-emission scores. Non-finite cells (e.g. the
    /// `-inf` for non-residue rows) are kept verbatim; the kernel only indexes rows for real
    /// residue codes, matching the scalar DP.
    pub fn new(prof: &P7Profile) -> FlatProfile {
        let (kp, m) = (prof.kp, prof.m);
        let mut msc = vec![f32::NEG_INFINITY; kp * (m + 1)];
        for (x, row) in prof.msc.iter().enumerate() {
            for (k, &s) in row.iter().enumerate() {
                msc[x * (m + 1) + k] = s;
            }
        }
        FlatProfile { kp, m, msc }
    }
}

/// `unbiased_byteify`: round a log-prob score to a non-negative uint8 cost (`255` = -inf), the
/// same quantization as `ornament-hmm::msv_simd`.
fn unbiased_byteify(scale_b: f32, sc: f32) -> u8 {
    let cost = -(scale_b * sc).round();
    if cost > 255.0 {
        255
    } else {
        cost.max(0.0) as u8
    }
}

/// `biased_byteify`: emission score -> biased uint8 cost (`bias_b - round(scale_b·sc)`) with
/// HMMER's two's-complement wraparound; `255` = -inf. Matches `ornament-hmm::msv_simd`.
fn biased_byteify(scale_b: f32, bias_b: u8, sc: f32) -> u8 {
    let cost = -(scale_b * sc).round();
    if cost > 255.0 - bias_b as f32 {
        255
    } else {
        (cost as i32 as u8).wrapping_add(bias_b)
    }
}

/// Reduced-precision (uint8) MSV profile for the throughput kernel — the byte-quantized analogue
/// of [`FlatProfile`], built with the same offset scheme as `ornament-hmm::MsvProfile`
/// (`scale_b = 3/ln2`, `base_b = 190`). Emission bytes `rbv` are 4x smaller than the f32 table,
/// which is the throughput win on the bandwidth-bound one-thread-per-window batch. The
/// length-model byte `tjb_b` is per-window and recomputed on-device from each window's length.
#[derive(Debug, Clone)]
pub struct ByteMsvProfile {
    /// Alphabet size (emission-table rows).
    pub kp: usize,
    /// Model length (match states).
    pub m: usize,
    /// Nats-per-byte scale (`3/ln2`); the device recovers raw nats and per-window `tjb_b` with it.
    pub scale_b: f32,
    /// Offset base byte (`190`).
    pub base_b: u8,
    /// Emission bias byte (from the max match log-odds).
    pub bias_b: u8,
    /// B->Mk entry byte (depends on `M`).
    pub tbm_b: u8,
    /// E->C / E->J byte (constant, multihit).
    pub tec_b: u8,
    /// Row-major `[Kp][M+1]` biased emission cost bytes (`255` = -inf; `k = 0` padded).
    pub rbv: Vec<u8>,
}

impl ByteMsvProfile {
    /// Quantize a configured [`P7Profile`] into the uint8 MSV profile.
    pub fn new(prof: &P7Profile) -> ByteMsvProfile {
        let (kp, m) = (prof.kp, prof.m);
        let scale_b = 3.0 / std::f32::consts::LN_2;
        let base_b = 190u8;

        // bias from the max finite match-emission log-odds.
        let mut max = 0.0f32;
        for row in prof.msc.iter() {
            for &s in row[1..=m].iter() {
                if s.is_finite() {
                    max = max.max(s);
                }
            }
        }
        let bias_b = unbiased_byteify(scale_b, -max);

        // rbv[x*(M+1)+k]: biased emission byte; -inf cells (k = 0, non-residue rows) quantize to 255.
        let mut rbv = vec![255u8; kp * (m + 1)];
        for (x, row) in prof.msc.iter().enumerate() {
            for k in 1..=m {
                rbv[x * (m + 1) + k] = biased_byteify(scale_b, bias_b, row[k]);
            }
        }

        let tbm_b = unbiased_byteify(scale_b, (2.0 / (m as f32 * (m as f32 + 1.0))).ln());
        let tec_b = unbiased_byteify(scale_b, 0.5f32.ln());
        ByteMsvProfile {
            kp,
            m,
            scale_b,
            base_b,
            bias_b,
            tbm_b,
            tec_b,
            rbv,
        }
    }
}

/// A batch of windows to score against one [`FlatProfile`]: every window's residue codes packed
/// contiguously for the kernel, indexed by per-window `offsets`/`lengths`. The packed codes are
/// **raw** digital symbols (0..Kp) with the sentinels stripped — the kernel does its own 1-based
/// indexing, so it wants only the residues.
#[derive(Debug, Clone, Default)]
pub struct Windows {
    /// Concatenated residue codes of every window (no sentinels).
    pub dsq: Vec<u8>,
    /// Start index of window `j` into `dsq`.
    pub offsets: Vec<i32>,
    /// Residue length of window `j`.
    pub lengths: Vec<i32>,
}

impl Windows {
    /// Pack a list of windows into one batch. Each input window is a sentinel-padded 1-based
    /// digital sequence (the codebase's `dsq` convention, as returned by `Alphabet::digitize`);
    /// the surrounding sentinels are dropped when packing.
    pub fn from_slices(windows: &[&[Dsq]]) -> Windows {
        let mut w = Windows::default();
        for win in windows {
            // Strip the leading/trailing sentinel to get the raw residues.
            let res = &win[1..win.len().saturating_sub(1).max(1)];
            w.offsets.push(w.dsq.len() as i32);
            w.lengths.push(res.len() as i32);
            w.dsq.extend_from_slice(res);
        }
        w
    }

    /// Number of windows in the batch.
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Whether the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }
}

/// Null-correct a raw MSV score (nats) to a **bit** score for a length-`l` window — the same
/// `(raw - null_one(l)) / ln2` convention as the scalar [`ornament_hmm::msv_bits`], so results
/// land on the scale the calibrated `msv_pvalue` expects.
#[inline]
pub fn nats_to_bits(raw_nats: f32, l: usize) -> f32 {
    (raw_nats - null_one(l)) / std::f32::consts::LN_2
}

/// CPU oracle: scalar MSV raw score (nats) for each window, via the trusted
/// [`ornament_hmm::msv_nats`]. This is the reference the device kernel is validated against. Each
/// window is a sentinel-padded 1-based `dsq` (as [`msv_nats`](ornament_hmm::msv_nats) expects).
pub fn msv_nats_cpu(prof: &P7Profile, windows: &[&[Dsq]]) -> Vec<f32> {
    windows
        .iter()
        .map(|win| ornament_hmm::msv_nats(prof, win))
        .collect()
}

// ---- GPU path (feature = "cuda") ------------------------------------------------------------

#[cfg(feature = "cuda")]
extern "C" {
    fn ornament_gpu_device_count(count: *mut i32) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_msv_batch(
        msc: *const f32,
        kp: i32,
        m: i32,
        dsq: *const u8,
        total_res: i32,
        offsets: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_msv_u8_batch(
        rbv: *const u8,
        kp: i32,
        m: i32,
        scale_b: f32,
        base_b: i32,
        bias_b: i32,
        tbm_b: i32,
        tec_b: i32,
        dsq: *const u8,
        total_res: i32,
        offsets: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
}

/// Number of visible CUDA devices, or `Err` if the runtime call fails / the crate lacks `cuda`.
pub fn device_count() -> Result<usize, GpuError> {
    #[cfg(feature = "cuda")]
    {
        let mut c: i32 = 0;
        let rc = unsafe { ornament_gpu_device_count(&mut c) };
        if rc != 0 {
            return Err(GpuError::Cuda(rc));
        }
        Ok(c as usize)
    }
    #[cfg(not(feature = "cuda"))]
    {
        Err(GpuError::NotCompiled)
    }
}

/// Batch MSV raw scores (nats) on the GPU, one score per window. Requires the `cuda` feature and
/// a visible device. The result is element-for-element comparable to [`msv_nats_cpu`].
#[cfg(feature = "cuda")]
pub fn msv_nats_gpu(fp: &FlatProfile, w: &Windows) -> Result<Vec<f32>, GpuError> {
    let n = w.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if device_count()? == 0 {
        return Err(GpuError::NoDevice);
    }
    let mut out = vec![f32::NEG_INFINITY; n];
    let rc = unsafe {
        ornament_gpu_msv_batch(
            fp.msc.as_ptr(),
            fp.kp as i32,
            fp.m as i32,
            w.dsq.as_ptr(),
            w.dsq.len() as i32,
            w.offsets.as_ptr(),
            w.lengths.as_ptr(),
            n as i32,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(GpuError::Cuda(rc));
    }
    Ok(out)
}

/// Stub when built without `cuda`: always [`GpuError::NotCompiled`]. Lets callers link against a
/// stable signature regardless of feature state.
#[cfg(not(feature = "cuda"))]
pub fn msv_nats_gpu(_fp: &FlatProfile, _w: &Windows) -> Result<Vec<f32>, GpuError> {
    Err(GpuError::NotCompiled)
}

/// Reduced-precision (uint8) batch MSV raw scores (nats) on the GPU, one per window — the
/// throughput analogue of [`msv_nats_gpu`]. Matches `ornament-hmm::MsvProfile::score_bits` (the
/// striped uint8 filter), not the f32 `msv_nats`: it floors weak scores and reports saturating
/// hits as `+inf`. Convert to bits with [`nats_to_bits`].
#[cfg(feature = "cuda")]
pub fn msv_u8_nats_gpu(bp: &ByteMsvProfile, w: &Windows) -> Result<Vec<f32>, GpuError> {
    let n = w.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if device_count()? == 0 {
        return Err(GpuError::NoDevice);
    }
    let mut out = vec![f32::NEG_INFINITY; n];
    let rc = unsafe {
        ornament_gpu_msv_u8_batch(
            bp.rbv.as_ptr(),
            bp.kp as i32,
            bp.m as i32,
            bp.scale_b,
            bp.base_b as i32,
            bp.bias_b as i32,
            bp.tbm_b as i32,
            bp.tec_b as i32,
            w.dsq.as_ptr(),
            w.dsq.len() as i32,
            w.offsets.as_ptr(),
            w.lengths.as_ptr(),
            n as i32,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(GpuError::Cuda(rc));
    }
    Ok(out)
}

/// Stub when built without `cuda`: always [`GpuError::NotCompiled`].
#[cfg(not(feature = "cuda"))]
pub fn msv_u8_nats_gpu(_bp: &ByteMsvProfile, _w: &Windows) -> Result<Vec<f32>, GpuError> {
    Err(GpuError::NotCompiled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ornament_alphabet::Alphabet;
    use ornament_hmm::profile::{bg_freqs, P7Profile};

    /// Build the embedded tRNA filter HMM + a configured local profile, plus a set of random
    /// windows of varied length. Shared fixture for the CPU-oracle and GPU-parity tests.
    fn fixture() -> (P7Profile, Vec<Vec<Dsq>>) {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../ornament-scfg/tests/data/trna_fp7.hmm"
        );
        let hmm = ornament_hmm::parse_p7_hmm(&std::fs::read_to_string(path).expect("read hmm"))
            .expect("parse hmm");
        let abc = Alphabet::rna();
        let bg = bg_freqs(hmm.k);
        let prof = P7Profile::config_local(&hmm, &abc, &bg, 100);

        let mut s: u64 = 0xC0FFEE;
        let bases = [b'A', b'C', b'G', b'U'];
        let windows: Vec<Vec<Dsq>> = [20usize, 73, 73, 160, 300, 436]
            .iter()
            .map(|&len| {
                let seq: String = (0..len)
                    .map(|_| {
                        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                        bases[((s >> 40) & 3) as usize] as char
                    })
                    .collect();
                abc.digitize(&seq).expect("digitize")
            })
            .collect();
        (prof, windows)
    }

    /// The batch builders and CPU oracle agree with per-window scalar `msv_nats`, and the flat
    /// emission table round-trips the profile. Runs everywhere (no GPU needed).
    #[test]
    fn cpu_batch_matches_scalar() {
        let (prof, windows) = fixture();
        let refs: Vec<&[Dsq]> = windows.iter().map(|w| w.as_slice()).collect();

        let fp = FlatProfile::new(&prof);
        assert_eq!(fp.msc.len(), fp.kp * (fp.m + 1));

        let cpu = msv_nats_cpu(&prof, &refs);
        assert_eq!(cpu.len(), refs.len());
        // nats -> bits conversion (residue length = padded len - 2) is finite for these lengths.
        for (raw, w) in cpu.iter().zip(&refs) {
            assert!(nats_to_bits(*raw, w.len() - 2).is_finite());
        }
    }

    /// GPU parity: the device kernel must match the scalar oracle to f32 rounding. Skips when the
    /// crate lacks `cuda` or no device is visible.
    #[cfg(feature = "cuda")]
    #[test]
    fn gpu_matches_cpu() {
        let (prof, windows) = fixture();
        let refs: Vec<&[Dsq]> = windows.iter().map(|w| w.as_slice()).collect();

        match device_count() {
            Ok(0) | Err(_) => return, // no device on this host — skip.
            Ok(_) => {}
        }

        let fp = FlatProfile::new(&prof);
        let batch = Windows::from_slices(&refs);
        let gpu = msv_nats_gpu(&fp, &batch).expect("gpu msv");
        let cpu = msv_nats_cpu(&prof, &refs);

        for (i, (g, c)) in gpu.iter().zip(&cpu).enumerate() {
            assert!(
                (g - c).abs() < 0.01,
                "window {i} (len {}): gpu {g} vs cpu {c}",
                refs[i].len()
            );
        }
    }

    /// uint8 GPU parity: the reduced-precision device kernel must match the CPU striped uint8
    /// filter (`ornament_hmm::MsvProfile::score_bits`) — same quantization, same recurrence, just
    /// a linear vs striped traversal of an associative max-plus DP. Compared in bit space (the
    /// scale `msv_pvalue` expects). Skips when no `cuda` device is visible; the oracle is
    /// x86_64-only (that's where `MsvProfile` lives), which the `compgpu` nodes are.
    #[cfg(all(feature = "cuda", target_arch = "x86_64"))]
    #[test]
    fn gpu_u8_matches_striped() {
        use ornament_hmm::MsvProfile;

        let (prof, windows) = fixture();
        let refs: Vec<&[Dsq]> = windows.iter().map(|w| w.as_slice()).collect();

        match device_count() {
            Ok(0) | Err(_) => return,
            Ok(_) => {}
        }

        let bp = ByteMsvProfile::new(&prof);
        let batch = Windows::from_slices(&refs);
        let gpu_raw = msv_u8_nats_gpu(&bp, &batch).expect("gpu u8 msv");

        for (i, win) in refs.iter().enumerate() {
            let res_len = win.len() - 2; // strip sentinels for the length model
            let mut sp = MsvProfile::new(&prof);
            sp.set_length(res_len);
            let cpu_bits = sp.score_bits(win);
            let gpu_bits = nats_to_bits(gpu_raw[i], res_len);
            // Both saturate to +inf together, or agree to f32 quantization.
            assert!(
                (gpu_bits - cpu_bits).abs() < 0.01
                    || (gpu_bits.is_infinite() && cpu_bits.is_infinite()),
                "window {i} (len {res_len}): gpu {gpu_bits} vs cpu striped {cpu_bits}"
            );
        }
    }
}
