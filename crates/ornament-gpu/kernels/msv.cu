// Batch MSV filter kernel + CUDA-Runtime-API host wrappers.
//
// This is the GPU port of the p7 MSV prefilter (ungapped multi-segment Viterbi), the
// cheapest-first stage of the CM acceleration cascade. It is a faithful transcription of the
// scalar max-plus recurrence in `ornament-hmm/src/msv.rs::msv_nats` — same semiring, same
// length model, same B->Mk entry — kept in **f32** (not the uint8 quantization) so the result
// matches the scalar oracle to f32 rounding and needs no calibration tuning.
//
// Parallelization: **one thread per window**. The all-model / many-genomes scan produces an
// enormous batch of independent (model, window) tiles; MSV has a sequential dependency down the
// sequence and a wavefront dependency across model positions, so intra-DP parallelism is
// awkward, but *across* windows it is embarrassingly parallel. One thread streams one window's
// DP end-to-end. A single model's emission table (msc) is shared read-only by every thread.
//
// Memory: each thread keeps one rolling DP row `mmx[0..=M]` (MSV needs only M(i-1,k-1), so a
// single row scanned k = M..1 suffices — mmx[k-1] still holds row i-1 when we write mmx[k]).
// That row lives in a caller-provided global scratch buffer (stride M+1 per thread) so M is
// unbounded; a production kernel would stage short models through shared memory.

#include <cstdint>
#include <math.h>

#ifndef NU
#define NU 2.0f          // expected hit count (HMMER MSV default)
#endif

#define NINF (-INFINITY)

__device__ __forceinline__ float fmaxf_(float a, float b) { return a > b ? a : b; }

// One thread scores window `j`. `msc` is row-major [Kp][M+1] (msc[x*(M+1)+k]); `dsq` is the
// concatenation of all windows' residue codes (0..Kp, no sentinels); window j occupies
// dsq[offsets[j] .. offsets[j]+lengths[j]). `scratch` gives each thread its (M+1) DP row.
// Output `out[j]` is the raw MSV score in nats (null correction / bit conversion is done host-
// side, exactly as the scalar `msv_bits` does).
extern "C" __global__ void msv_batch_kernel(
    const float*   __restrict__ msc,
    int                          Kp,
    int                          M,
    const uint8_t* __restrict__ dsq,
    const int*     __restrict__ offsets,
    const int*     __restrict__ lengths,
    int                          n,
    float*         __restrict__ scratch,
    float*         __restrict__ out)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;

    int   L    = lengths[j];
    int   base = offsets[j];
    float* mmx = scratch + (size_t)j * (M + 1);

    // MSV-specific transitions (generic_msv.c). tbmk/tej/tec depend only on M; tloop/tmove on L.
    float tloop = logf((float)L / ((float)L + 3.0f));
    float tmove = logf(3.0f / ((float)L + 3.0f));
    float tbmk  = logf(2.0f / ((float)M * ((float)M + 1.0f)));
    float tej   = logf((NU - 1.0f) / NU);
    float tec   = logf(1.0f / NU);

    for (int k = 0; k <= M; ++k) mmx[k] = NINF;

    float xN = 0.0f, xB = tmove, xJ = NINF, xC = NINF;

    for (int i = 1; i <= L; ++i) {
        int          x   = (int)dsq[base + (i - 1)];
        const float* row = msc + (size_t)x * (M + 1);
        float        xE  = NINF;
        float        xBt = xB + tbmk;   // shared B->Mk entry for this row

        // k = M..1 so mmx[k-1] is still the previous row's value when we overwrite mmx[k].
        for (int k = M; k >= 1; --k) {
            float mk = row[k] + fmaxf_(mmx[k - 1], xBt);
            mmx[k] = mk;
            xE = fmaxf_(xE, mk);
        }

        xJ = fmaxf_(xJ + tloop, xE + tej);
        xC = fmaxf_(xC + tloop, xE + tec);
        xN = xN + tloop;
        xB = fmaxf_(xN + tmove, xJ + tmove);
    }

    out[j] = xC + tmove;   // L == 0 -> xC stays NINF -> -inf, matching the scalar guard
}

// ---- Reduced-precision (uint8) MSV ----------------------------------------------------------
// Throughput variant: the same recurrence in HMMER's reduced-precision offset-uint8 arithmetic
// (impl_sse/msvfilter.c, mirrored by `ornament-hmm::msv_simd`). Scores live in offset "third-bit"
// units; a strong hit *saturates* the uint8 range and is reported as +inf so a higher-precision
// filter re-scores it. The emission table is 4x smaller than f32 (better bandwidth/cache), which
// is the whole point on the memory-bound one-thread-per-window batch.
//
// One thread per window, linear scan k = M..1: because MSV is an associative max-plus recurrence,
// the linear order computes the *same* bytes as the CPU's striped kernel, so this matches
// `MsvProfile::score_bits` (not the f32 `msv_nats`). Saturating ops are done in `int` and clamped.

__device__ __forceinline__ int sadd_u8(int a, int b) { int s = a + b; return s > 255 ? 255 : s; }
__device__ __forceinline__ int ssub_u8(int a, int b) { int s = a - b; return s < 0 ? 0 : s; }

// tjb_b = unbiased_byteify(scale_b, log(3/(L+3))): round(-scale*sc), clamped to [0,255].
__device__ __forceinline__ int tjb_byte(float scale_b, int L) {
    float sc   = logf(3.0f / ((float)L + 3.0f));
    float cost = -(scale_b * sc);
    cost = rintf(cost);
    if (cost > 255.0f) return 255;
    if (cost < 0.0f)   return 0;
    return (int)cost;
}

// rbv: biased emission cost bytes, row-major [Kp][M+1] (rbv[x*(M+1)+k]); byte 255 = -inf cost, and
// k = 0 is padded 255. base_b/bias_b/tbm_b/tec_b are the profile's constant bytes; scale_b recovers
// nats. Output is the raw MSV score in nats (host null-corrects to bits), or +inf on saturation.
extern "C" __global__ void msv_u8_batch_kernel(
    const uint8_t* __restrict__ rbv,
    int                          Kp,
    int                          M,
    float                        scale_b,
    int                          base_b,
    int                          bias_b,
    int                          tbm_b,
    int                          tec_b,
    const uint8_t* __restrict__ dsq,
    const int*     __restrict__ offsets,
    const int*     __restrict__ lengths,
    int                          n,
    uint8_t*       __restrict__ scratch,
    float*         __restrict__ out)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;

    int      L    = lengths[j];
    int      base = offsets[j];
    uint8_t* mmx  = scratch + (size_t)j * (M + 1);

    int tjb  = tjb_byte(scale_b, L);
    int tjbm = (tjb + tbm_b) & 0xff;   // HMMER wraps this byte sum, then saturating-subtracts it

    for (int k = 0; k <= M; ++k) mmx[k] = 0;   // 0 == the -inf-equivalent floor in the offset scheme

    int xJ = 0;
    int xB = ssub_u8(base_b, tjbm);

    for (int i = 1; i <= L; ++i) {
        int            x   = (int)dsq[base + (i - 1)];
        const uint8_t* row = rbv + (size_t)x * (M + 1);
        int            xE  = 0;

        for (int k = M; k >= 1; --k) {
            int sv = mmx[k - 1] > xB ? mmx[k - 1] : xB;   // max(M(i-1,k-1), xB)
            sv = sadd_u8(sv, bias_b);
            sv = ssub_u8(sv, (int)row[k]);
            mmx[k] = (uint8_t)sv;
            if (sv > xE) xE = sv;
        }

        if (sadd_u8(xE, bias_b) == 255) { out[j] = INFINITY; return; }   // saturation -> pass

        int xJn = ssub_u8(xE, tec_b);
        if (xJn > xJ) xJ = xJn;
        int bmax = base_b > xJ ? base_b : xJ;
        xB = ssub_u8(bmax, tjbm);
    }

    // Recover raw nats (HMMER MSVFilter recovery); host converts to bits with null_one(L).
    out[j] = (((float)xJ - (float)tjb) - (float)base_b) / scale_b - 3.0f;
}

// ---- Host wrappers (CUDA Runtime API) -------------------------------------------------------
// Thin, synchronous launch helpers callable over FFI. They own the device allocations for one
// batch call: upload inputs, launch, copy scores back, free. Return 0 on success, or the
// cudaError_t (>0) of the first failing call so the Rust side can surface it.

static int check(cudaError_t e) { return (int)e; }

#define CUDA_TRY(call) do { cudaError_t _e = (call); if (_e != cudaSuccess) { rc = (int)_e; goto cleanup; } } while (0)

extern "C" int ornament_gpu_device_count(int* count) {
    return check(cudaGetDeviceCount(count));
}

// ---- Resident strand ------------------------------------------------------------------------
// Upload a strand's residues to the device ONCE and score many windows against it as (start,len)
// offsets — instead of re-copying overlapping window bytes per batch. The scan tiles at length 2W
// step W (2x overlap), so packing windows moves each base ~2x over PCIe; a resident strand moves
// it once, and the same device buffer is reused across the MSV/Viterbi/Forward cascade.

// Upload `len` residue bytes; return the device pointer in `*d_out` (opaque handle for Rust).
extern "C" int ornament_gpu_strand_upload(const uint8_t* dsq, int len, void** d_out) {
    int rc = 0;
    uint8_t* d = nullptr;
    *d_out = nullptr;
    if (len < 0) return (int)cudaErrorInvalidValue;
    CUDA_TRY(cudaMalloc(&d, (size_t)len * sizeof(uint8_t)));
    CUDA_TRY(cudaMemcpy(d, dsq, (size_t)len * sizeof(uint8_t), cudaMemcpyHostToDevice));
    *d_out = d;
    return 0;
cleanup:
    if (d) cudaFree(d);
    return rc;
}

// Free a device strand from `ornament_gpu_strand_upload`.
extern "C" int ornament_gpu_strand_free(void* d) {
    return check(cudaFree(d));
}

// f32 MSV over a resident strand. `d_strand` is a device pointer (from strand_upload); `starts`
// are 0-based offsets into it and `lengths` the window lengths. msc: Kp*(M+1) floats (host).
extern "C" int ornament_gpu_msv_batch_resident(
    const float*   msc,
    int            Kp,
    int            M,
    const void*    d_strand,
    const int*     starts,
    const int*     lengths,
    int            n,
    float*         out)
{
    int rc = 0;
    float *d_msc = nullptr, *d_scratch = nullptr, *d_out = nullptr;
    int   *d_off = nullptr, *d_len = nullptr;

    if (n <= 0) return 0;

    size_t msc_bytes     = (size_t)Kp * (M + 1) * sizeof(float);
    size_t scratch_bytes = (size_t)n * (M + 1) * sizeof(float);

    CUDA_TRY(cudaMalloc(&d_msc,     msc_bytes));
    CUDA_TRY(cudaMalloc(&d_off,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_len,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_out,     (size_t)n * sizeof(float)));
    CUDA_TRY(cudaMalloc(&d_scratch, scratch_bytes));

    CUDA_TRY(cudaMemcpy(d_msc, msc, msc_bytes,                   cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_off, starts,  (size_t)n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_len, lengths, (size_t)n * sizeof(int), cudaMemcpyHostToDevice));

    {
        int threads = 128;
        int blocks  = (n + threads - 1) / threads;
        msv_batch_kernel<<<blocks, threads>>>(
            d_msc, Kp, M, (const uint8_t*)d_strand, d_off, d_len, n, d_scratch, d_out);
        CUDA_TRY(cudaGetLastError());
        CUDA_TRY(cudaDeviceSynchronize());
    }

    CUDA_TRY(cudaMemcpy(out, d_out, (size_t)n * sizeof(float), cudaMemcpyDeviceToHost));

cleanup:
    cudaFree(d_msc);
    cudaFree(d_off);
    cudaFree(d_len);
    cudaFree(d_out);
    cudaFree(d_scratch);
    return rc;
}

// Reduced-precision (uint8) MSV over a resident strand. rbv: Kp*(M+1) bytes (host); scale/base/
// bias/tbm/tec from the host `ByteMsvProfile`. out: n floats (raw nats / +inf).
extern "C" int ornament_gpu_msv_u8_batch_resident(
    const uint8_t* rbv,
    int            Kp,
    int            M,
    float          scale_b,
    int            base_b,
    int            bias_b,
    int            tbm_b,
    int            tec_b,
    const void*    d_strand,
    const int*     starts,
    const int*     lengths,
    int            n,
    float*         out)
{
    int rc = 0;
    uint8_t *d_rbv = nullptr, *d_scratch = nullptr;
    int     *d_off = nullptr, *d_len = nullptr;
    float   *d_out = nullptr;

    if (n <= 0) return 0;

    size_t rbv_bytes     = (size_t)Kp * (M + 1) * sizeof(uint8_t);
    size_t scratch_bytes = (size_t)n * (M + 1) * sizeof(uint8_t);

    CUDA_TRY(cudaMalloc(&d_rbv,     rbv_bytes));
    CUDA_TRY(cudaMalloc(&d_off,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_len,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_out,     (size_t)n * sizeof(float)));
    CUDA_TRY(cudaMalloc(&d_scratch, scratch_bytes));

    CUDA_TRY(cudaMemcpy(d_rbv, rbv, rbv_bytes,                   cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_off, starts,  (size_t)n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_len, lengths, (size_t)n * sizeof(int), cudaMemcpyHostToDevice));

    {
        int threads = 128;
        int blocks  = (n + threads - 1) / threads;
        msv_u8_batch_kernel<<<blocks, threads>>>(
            d_rbv, Kp, M, scale_b, base_b, bias_b, tbm_b, tec_b,
            (const uint8_t*)d_strand, d_off, d_len, n, d_scratch, d_out);
        CUDA_TRY(cudaGetLastError());
        CUDA_TRY(cudaDeviceSynchronize());
    }

    CUDA_TRY(cudaMemcpy(out, d_out, (size_t)n * sizeof(float), cudaMemcpyDeviceToHost));

cleanup:
    cudaFree(d_rbv);
    cudaFree(d_off);
    cudaFree(d_len);
    cudaFree(d_out);
    cudaFree(d_scratch);
    return rc;
}
