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

// ---- Host wrappers (CUDA Runtime API) -------------------------------------------------------
// Thin, synchronous launch helpers callable over FFI. They own the device allocations for one
// batch call: upload inputs, launch, copy scores back, free. Return 0 on success, or the
// cudaError_t (>0) of the first failing call so the Rust side can surface it.

static int check(cudaError_t e) { return (int)e; }

#define CUDA_TRY(call) do { cudaError_t _e = (call); if (_e != cudaSuccess) { rc = (int)_e; goto cleanup; } } while (0)

extern "C" int ornament_gpu_device_count(int* count) {
    return check(cudaGetDeviceCount(count));
}

// msc: Kp*(M+1) floats. dsq: total_res u8. offsets/lengths: n ints. out: n floats (host).
extern "C" int ornament_gpu_msv_batch(
    const float*   msc,
    int            Kp,
    int            M,
    const uint8_t* dsq,
    int            total_res,
    const int*     offsets,
    const int*     lengths,
    int            n,
    float*         out)
{
    int rc = 0;
    float   *d_msc = nullptr, *d_scratch = nullptr, *d_out = nullptr;
    uint8_t *d_dsq = nullptr;
    int     *d_off = nullptr, *d_len = nullptr;

    if (n <= 0) return 0;

    size_t msc_bytes     = (size_t)Kp * (M + 1) * sizeof(float);
    size_t scratch_bytes = (size_t)n * (M + 1) * sizeof(float);

    CUDA_TRY(cudaMalloc(&d_msc,     msc_bytes));
    CUDA_TRY(cudaMalloc(&d_dsq,     (size_t)total_res * sizeof(uint8_t)));
    CUDA_TRY(cudaMalloc(&d_off,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_len,     (size_t)n * sizeof(int)));
    CUDA_TRY(cudaMalloc(&d_out,     (size_t)n * sizeof(float)));
    CUDA_TRY(cudaMalloc(&d_scratch, scratch_bytes));

    CUDA_TRY(cudaMemcpy(d_msc, msc, msc_bytes,                      cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_dsq, dsq, (size_t)total_res,              cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_off, offsets, (size_t)n * sizeof(int),    cudaMemcpyHostToDevice));
    CUDA_TRY(cudaMemcpy(d_len, lengths, (size_t)n * sizeof(int),    cudaMemcpyHostToDevice));

    {
        int threads = 128;
        int blocks  = (n + threads - 1) / threads;
        msv_batch_kernel<<<blocks, threads>>>(d_msc, Kp, M, d_dsq, d_off, d_len, n, d_scratch, d_out);
        CUDA_TRY(cudaGetLastError());
        CUDA_TRY(cudaDeviceSynchronize());
    }

    CUDA_TRY(cudaMemcpy(out, d_out, (size_t)n * sizeof(float), cudaMemcpyDeviceToHost));

cleanup:
    cudaFree(d_msc);
    cudaFree(d_dsq);
    cudaFree(d_off);
    cudaFree(d_len);
    cudaFree(d_out);
    cudaFree(d_scratch);
    return rc;
}
