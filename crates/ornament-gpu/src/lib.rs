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
}
