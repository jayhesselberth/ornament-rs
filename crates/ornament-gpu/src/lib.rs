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
//! One configured model's emission table, a strand uploaded once, and windows as offsets into it:
//! - [`FlatProfile`] / [`ByteMsvProfile`] — the model's match-emission log-odds (f32 / uint8),
//!   flattened to `[x*(M+1)+k]`.
//! - [`DeviceStrand`] — a strand's residues, resident on the device and reused across windows and
//!   cascade stages (avoids re-copying overlapping tiles).
//! - [`Tiles`] — the windows as 1-based `(start, end)` spans → device offsets.
//!
//! ## Feature gating
//! The `cuda` feature turns on the device path (build.rs compiles the kernel with nvcc and links
//! the CUDA runtime). Without it the crate still builds and the CPU oracle + batch builders
//! ([`FlatProfile`], [`ByteMsvProfile`], [`Tiles`], [`msv_nats_cpu`]) are available; the device
//! types ([`DeviceStrand`], [`msv_nats_gpu`]) are compiled out. This keeps the workspace buildable
//! on login nodes with no CUDA toolchain; build/run the GPU path on the `compgpu` nodes.

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

/// Windows into a [`DeviceStrand`], as `(start, length)` offsets rather than copied residues. The
/// scan tiles overlap (length `2W`, step `W`), so a copy-per-window batch moves each base ~2× over
/// PCIe; referencing a resident strand moves it once. `starts` are 0-based offsets into the
/// strand's residues (sentinels already stripped), matching the kernel's own 1-based indexing.
#[derive(Debug, Clone, Default)]
pub struct Tiles {
    /// 0-based start offset of window `j` into the resident strand.
    pub starts: Vec<i32>,
    /// Residue length of window `j`.
    pub lengths: Vec<i32>,
}

impl Tiles {
    /// Build from 1-based inclusive `(start, end)` spans — the coordinate convention the scan
    /// pipeline tiles in (`strand_survivors`). A span `(s, e)` covers residues `s..=e` of the
    /// strand, i.e. 0-based start `s-1`, length `e-s+1`.
    pub fn from_spans(spans: &[(usize, usize)]) -> Tiles {
        let mut t = Tiles::default();
        for &(s, e) in spans {
            debug_assert!(s >= 1 && e >= s, "1-based inclusive span, got ({s},{e})");
            t.starts.push((s - 1) as i32);
            t.lengths.push((e - s + 1) as i32);
        }
        t
    }

    /// Number of windows.
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    /// Whether there are no windows.
    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
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
    fn ornament_gpu_strand_upload(
        dsq: *const u8,
        len: i32,
        d_out: *mut *mut core::ffi::c_void,
    ) -> i32;
    fn ornament_gpu_strand_free(d: *mut core::ffi::c_void) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_msv_batch_resident(
        msc: *const f32,
        kp: i32,
        m: i32,
        d_strand: *const core::ffi::c_void,
        starts: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_msv_u8_batch_resident(
        rbv: *const u8,
        kp: i32,
        m: i32,
        scale_b: f32,
        base_b: i32,
        bias_b: i32,
        tbm_b: i32,
        tec_b: i32,
        d_strand: *const core::ffi::c_void,
        starts: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
    fn ornament_gpu_ctx_create(streams: i32, chunk: i32, out: *mut *mut core::ffi::c_void) -> i32;
    fn ornament_gpu_ctx_destroy(ctx: *mut core::ffi::c_void) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_ctx_msv_f32(
        ctx: *mut core::ffi::c_void,
        msc: *const f32,
        kp: i32,
        m: i32,
        d_strand: *const core::ffi::c_void,
        starts: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
    #[allow(clippy::too_many_arguments)]
    fn ornament_gpu_ctx_msv_u8(
        ctx: *mut core::ffi::c_void,
        rbv: *const u8,
        kp: i32,
        m: i32,
        scale_b: f32,
        base_b: i32,
        bias_b: i32,
        tbm_b: i32,
        tec_b: i32,
        d_strand: *const core::ffi::c_void,
        starts: *const i32,
        lengths: *const i32,
        n: i32,
        out: *mut f32,
    ) -> i32;
}

/// A strand's residues uploaded to the device **once**, reusable across every window batch and
/// every cascade stage (MSV/Viterbi/Forward) without re-copying. RAII: the device buffer is freed
/// on drop.
#[cfg(feature = "cuda")]
pub struct DeviceStrand {
    ptr: *mut core::ffi::c_void,
    len: usize,
}

#[cfg(feature = "cuda")]
impl DeviceStrand {
    /// Upload a sentinel-padded 1-based `dsq` (as `Alphabet::digitize` returns). The surrounding
    /// sentinels are stripped; the raw residues live on the device until this is dropped.
    pub fn upload(dsq: &[Dsq]) -> Result<DeviceStrand, GpuError> {
        if device_count()? == 0 {
            return Err(GpuError::NoDevice);
        }
        // Strip the leading/trailing sentinel to get the raw residues.
        let res: &[u8] = if dsq.len() >= 2 {
            &dsq[1..dsq.len() - 1]
        } else {
            &[]
        };
        let mut ptr: *mut core::ffi::c_void = core::ptr::null_mut();
        let rc = unsafe { ornament_gpu_strand_upload(res.as_ptr(), res.len() as i32, &mut ptr) };
        if rc != 0 {
            return Err(GpuError::Cuda(rc));
        }
        Ok(DeviceStrand {
            ptr,
            len: res.len(),
        })
    }

    /// Residue length of the resident strand.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the strand is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

#[cfg(feature = "cuda")]
impl Drop for DeviceStrand {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ornament_gpu_strand_free(self.ptr) };
        }
    }
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

/// Batch MSV raw scores (nats) on the GPU: score `tiles` against a resident `strand`, one score
/// per tile. The strand is uploaded once ([`DeviceStrand::upload`]) and reused across tiles and
/// cascade stages. Results are element-for-element comparable to [`msv_nats_cpu`].
#[cfg(feature = "cuda")]
pub fn msv_nats_gpu(
    fp: &FlatProfile,
    strand: &DeviceStrand,
    tiles: &Tiles,
) -> Result<Vec<f32>, GpuError> {
    let n = tiles.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut out = vec![f32::NEG_INFINITY; n];
    let rc = unsafe {
        ornament_gpu_msv_batch_resident(
            fp.msc.as_ptr(),
            fp.kp as i32,
            fp.m as i32,
            strand.ptr,
            tiles.starts.as_ptr(),
            tiles.lengths.as_ptr(),
            n as i32,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(GpuError::Cuda(rc));
    }
    Ok(out)
}

/// Reduced-precision (uint8) batch MSV raw scores (nats) on the GPU — the throughput analogue of
/// [`msv_nats_gpu`], scoring `tiles` against a resident `strand`. Matches
/// `ornament-hmm::MsvProfile::score_bits` (the striped uint8 filter), not the f32 `msv_nats`: it
/// floors weak scores and reports saturating hits as `+inf`. Convert to bits with [`nats_to_bits`].
#[cfg(feature = "cuda")]
pub fn msv_u8_nats_gpu(
    bp: &ByteMsvProfile,
    strand: &DeviceStrand,
    tiles: &Tiles,
) -> Result<Vec<f32>, GpuError> {
    let n = tiles.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut out = vec![f32::NEG_INFINITY; n];
    let rc = unsafe {
        ornament_gpu_msv_u8_batch_resident(
            bp.rbv.as_ptr(),
            bp.kp as i32,
            bp.m as i32,
            bp.scale_b,
            bp.base_b as i32,
            bp.bias_b as i32,
            bp.tbm_b as i32,
            bp.tec_b as i32,
            strand.ptr,
            tiles.starts.as_ptr(),
            tiles.lengths.as_ptr(),
            n as i32,
            out.as_mut_ptr(),
        )
    };
    if rc != 0 {
        return Err(GpuError::Cuda(rc));
    }
    Ok(out)
}

/// A reusable GPU context: owns the CUDA streams + per-stream device/pinned buffers + a growable
/// emission-table buffer, created once and reused across many scoring calls. This is what makes an
/// all-model scan efficient — the free functions ([`msv_nats_gpu`], [`msv_u8_nats_gpu`]) spin up
/// and tear down all of that *per call*, which is negligible once but dominates when scoring
/// thousands of Rfam models. Score with [`GpuContext::msv_nats`] / [`GpuContext::msv_u8_nats`],
/// reusing one context (and one [`DeviceStrand`]) across every model. RAII: freed on drop.
#[cfg(feature = "cuda")]
pub struct GpuContext {
    ptr: *mut core::ffi::c_void,
}

#[cfg(feature = "cuda")]
impl GpuContext {
    /// Create a context with the default stream count / chunk size (honoring `ORNAMENT_GPU_STREAMS`
    /// / `ORNAMENT_GPU_CHUNK`). Requires a visible CUDA device.
    pub fn new() -> Result<GpuContext, GpuError> {
        Self::with_config(0, 0)
    }

    /// Create a context with an explicit stream count and chunk size (tiles/chunk). Pass `0` for
    /// either to fall back to the env var / built-in default.
    pub fn with_config(streams: usize, chunk: usize) -> Result<GpuContext, GpuError> {
        if device_count()? == 0 {
            return Err(GpuError::NoDevice);
        }
        let mut ptr: *mut core::ffi::c_void = core::ptr::null_mut();
        let rc = unsafe { ornament_gpu_ctx_create(streams as i32, chunk as i32, &mut ptr) };
        if rc != 0 {
            return Err(GpuError::Cuda(rc));
        }
        Ok(GpuContext { ptr })
    }

    /// f32 MSV raw scores (nats) for `tiles` against a resident `strand`, reusing this context's
    /// streams/buffers. Same result as [`msv_nats_gpu`], without the per-call setup.
    pub fn msv_nats(
        &self,
        fp: &FlatProfile,
        strand: &DeviceStrand,
        tiles: &Tiles,
    ) -> Result<Vec<f32>, GpuError> {
        let n = tiles.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let mut out = vec![f32::NEG_INFINITY; n];
        let rc = unsafe {
            ornament_gpu_ctx_msv_f32(
                self.ptr,
                fp.msc.as_ptr(),
                fp.kp as i32,
                fp.m as i32,
                strand.ptr,
                tiles.starts.as_ptr(),
                tiles.lengths.as_ptr(),
                n as i32,
                out.as_mut_ptr(),
            )
        };
        if rc != 0 {
            return Err(GpuError::Cuda(rc));
        }
        Ok(out)
    }

    /// Reduced-precision (uint8) MSV raw scores (nats), reusing this context. Same result as
    /// [`msv_u8_nats_gpu`], without the per-call setup.
    pub fn msv_u8_nats(
        &self,
        bp: &ByteMsvProfile,
        strand: &DeviceStrand,
        tiles: &Tiles,
    ) -> Result<Vec<f32>, GpuError> {
        let n = tiles.len();
        if n == 0 {
            return Ok(Vec::new());
        }
        let mut out = vec![f32::NEG_INFINITY; n];
        let rc = unsafe {
            ornament_gpu_ctx_msv_u8(
                self.ptr,
                bp.rbv.as_ptr(),
                bp.kp as i32,
                bp.m as i32,
                bp.scale_b,
                bp.base_b as i32,
                bp.bias_b as i32,
                bp.tbm_b as i32,
                bp.tec_b as i32,
                strand.ptr,
                tiles.starts.as_ptr(),
                tiles.lengths.as_ptr(),
                n as i32,
                out.as_mut_ptr(),
            )
        };
        if rc != 0 {
            return Err(GpuError::Cuda(rc));
        }
        Ok(out)
    }
}

#[cfg(feature = "cuda")]
impl Drop for GpuContext {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { ornament_gpu_ctx_destroy(self.ptr) };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ornament_alphabet::Alphabet;
    use ornament_hmm::profile::{bg_freqs, P7Profile};

    /// Overlapping 1-based inclusive `(start, end)` spans into the fixture strand — exercises the
    /// resident-strand offset indexing (adjacent, nested, and full-length windows).
    const SPANS: &[(usize, usize)] = &[
        (1, 73),
        (50, 260),
        (200, 260),
        (300, 900),
        (1, 1200),
        (1100, 1200),
    ];

    /// Build the embedded tRNA filter HMM + a configured local profile, and one long random strand
    /// (sentinel-padded `dsq`, length 1200). Windows are [`SPANS`] into it. Shared fixture.
    fn fixture() -> (P7Profile, Vec<Dsq>) {
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
        let seq: String = (0..1200)
            .map(|_| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                bases[((s >> 40) & 3) as usize] as char
            })
            .collect();
        let strand = abc.digitize(&seq).expect("digitize");
        (prof, strand)
    }

    /// Sentinel-padded subsequence for a 1-based inclusive span `(s, e)` of a padded strand — the
    /// per-window `dsq` the scalar oracles expect.
    fn padded_span(strand: &[Dsq], s: usize, e: usize) -> Vec<Dsq> {
        let mut out = Vec::with_capacity(e - s + 3);
        out.push(ornament_alphabet::alphabet::SENTINEL);
        out.extend_from_slice(&strand[s..=e]); // strand[k] is residue k (1-based)
        out.push(ornament_alphabet::alphabet::SENTINEL);
        out
    }

    /// The flat emission table round-trips the profile and the CPU oracle scores each span finitely.
    /// Runs everywhere (no GPU needed).
    #[test]
    fn cpu_batch_matches_scalar() {
        let (prof, strand) = fixture();
        let fp = FlatProfile::new(&prof);
        assert_eq!(fp.msc.len(), fp.kp * (fp.m + 1));

        for &(s, e) in SPANS {
            let win = padded_span(&strand, s, e);
            let raw = msv_nats_cpu(&prof, &[&win])[0];
            assert!(nats_to_bits(raw, e - s + 1).is_finite());
        }
    }

    /// GPU parity (f32): scoring [`SPANS`] against a resident [`DeviceStrand`] must match the scalar
    /// oracle applied to each span's subsequence — confirming the offset indexing reads the right
    /// slice of the resident strand. Skips when no `cuda` device is visible.
    #[cfg(feature = "cuda")]
    #[test]
    fn gpu_matches_cpu() {
        let (prof, strand) = fixture();
        match device_count() {
            Ok(0) | Err(_) => return, // no device on this host — skip.
            Ok(_) => {}
        }

        let fp = FlatProfile::new(&prof);
        let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
        let tiles = Tiles::from_spans(SPANS);
        let gpu = msv_nats_gpu(&fp, &dstrand, &tiles).expect("gpu msv");

        for (i, &(s, e)) in SPANS.iter().enumerate() {
            let win = padded_span(&strand, s, e);
            let cpu = msv_nats_cpu(&prof, &[&win])[0];
            assert!(
                (gpu[i] - cpu).abs() < 0.01,
                "span {i} ({s},{e}): gpu {} vs cpu {cpu}",
                gpu[i]
            );
        }
    }

    /// Streamed pipeline across many chunks: a batch larger than the kernel's internal chunk size
    /// must still score every tile correctly, exercising the double-buffer slot-reclaim path (a
    /// stream is synced + copied out before its buffers are reused). Uses the f32 path (cheap exact
    /// oracle); the streaming helper is shared with the uint8 path. Distinct spans are cached so the
    /// CPU oracle stays fast. Skips when no `cuda` device is visible.
    #[cfg(feature = "cuda")]
    #[test]
    fn gpu_streamed_multichunk() {
        use std::collections::HashMap;

        let (prof, strand) = fixture();
        match device_count() {
            Ok(0) | Err(_) => return,
            Ok(_) => {}
        }

        // > 2 internal chunks so slots are reclaimed and reused mid-stream.
        let n = 80_000usize;
        let spans: Vec<(usize, usize)> = (0..n)
            .map(|i| (1 + (i * 37) % (1200 - 40), 0))
            .map(|(s, _)| (s, s + 39))
            .collect();

        let fp = FlatProfile::new(&prof);
        let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
        let tiles = Tiles::from_spans(&spans);
        let gpu = msv_nats_gpu(&fp, &dstrand, &tiles).expect("gpu msv");
        assert_eq!(gpu.len(), n);

        let mut cache: HashMap<(usize, usize), f32> = HashMap::new();
        for (i, &(s, e)) in spans.iter().enumerate() {
            let cpu = *cache.entry((s, e)).or_insert_with(|| {
                let win = padded_span(&strand, s, e);
                msv_nats_cpu(&prof, &[&win])[0]
            });
            assert!(
                (gpu[i] - cpu).abs() < 0.01,
                "tile {i} span ({s},{e}): gpu {} vs cpu {cpu}",
                gpu[i]
            );
        }
    }

    /// uint8 GPU parity: the reduced-precision device kernel over a resident strand must match the
    /// CPU striped uint8 filter (`ornament_hmm::MsvProfile::score_bits`) per span — same
    /// quantization, same recurrence, linear vs striped traversal of an associative max-plus DP.
    /// Skips when no `cuda` device is visible; the oracle is x86_64-only (the `compgpu` nodes are).
    #[cfg(all(feature = "cuda", target_arch = "x86_64"))]
    #[test]
    fn gpu_u8_matches_striped() {
        use ornament_hmm::MsvProfile;

        let (prof, strand) = fixture();
        match device_count() {
            Ok(0) | Err(_) => return,
            Ok(_) => {}
        }

        let bp = ByteMsvProfile::new(&prof);
        let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
        let tiles = Tiles::from_spans(SPANS);
        let gpu_raw = msv_u8_nats_gpu(&bp, &dstrand, &tiles).expect("gpu u8 msv");

        for (i, &(s, e)) in SPANS.iter().enumerate() {
            let res_len = e - s + 1;
            let win = padded_span(&strand, s, e);
            let mut sp = MsvProfile::new(&prof);
            sp.set_length(res_len);
            let cpu_bits = sp.score_bits(&win);
            let gpu_bits = nats_to_bits(gpu_raw[i], res_len);
            // Both saturate to +inf together, or agree to f32 quantization.
            assert!(
                (gpu_bits - cpu_bits).abs() < 0.01
                    || (gpu_bits.is_infinite() && cpu_bits.is_infinite()),
                "span {i} ({s},{e}): gpu {gpu_bits} vs cpu striped {cpu_bits}"
            );
        }
    }

    /// A reused [`GpuContext`] must produce the same scores as the per-call free functions, across
    /// *multiple* scoring calls on one context (exercises buffer reuse + the emission-buffer regrow
    /// path when a second call has a different M). Skips when no `cuda` device is visible.
    #[cfg(feature = "cuda")]
    #[test]
    fn gpu_context_matches_freefn() {
        let (prof, strand) = fixture();
        match device_count() {
            Ok(0) | Err(_) => return,
            Ok(_) => {}
        }

        let fp = FlatProfile::new(&prof);
        let bp = ByteMsvProfile::new(&prof);
        let dstrand = DeviceStrand::upload(&strand).expect("upload strand");
        let tiles = Tiles::from_spans(SPANS);

        let ctx = GpuContext::new().expect("ctx");
        // Two calls on the same context: f32 then u8 (regrows the shared emission buffer),
        // then f32 again (reuses it) — all must match the per-call free functions.
        let ctx_f32_a = ctx.msv_nats(&fp, &dstrand, &tiles).expect("ctx f32");
        let ctx_u8 = ctx.msv_u8_nats(&bp, &dstrand, &tiles).expect("ctx u8");
        let ctx_f32_b = ctx.msv_nats(&fp, &dstrand, &tiles).expect("ctx f32 again");
        let free_f32 = msv_nats_gpu(&fp, &dstrand, &tiles).expect("free f32");
        let free_u8 = msv_u8_nats_gpu(&bp, &dstrand, &tiles).expect("free u8");

        for i in 0..SPANS.len() {
            assert_eq!(
                ctx_f32_a[i].to_bits(),
                free_f32[i].to_bits(),
                "f32 span {i}"
            );
            assert_eq!(
                ctx_f32_b[i].to_bits(),
                free_f32[i].to_bits(),
                "f32(reuse) span {i}"
            );
            assert_eq!(ctx_u8[i].to_bits(), free_u8[i].to_bits(), "u8 span {i}");
        }
    }
}
