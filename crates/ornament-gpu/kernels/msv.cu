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
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

// Shared-memory uint8 MSV. Identical recurrence to `msv_u8_batch_kernel`, but the per-thread DP
// row and the (small) emission table live in **shared** memory instead of global. That's the whole
// game: the global-scratch kernel is latency-bound because it hits the DP row M*L times over a
// serial dependency at ~400-cycle global latency; shared memory is ~30 cycles, so the same DP runs
// far faster. Requires (Kp + blockDim)*(M+1) bytes of dynamic shared memory — the host launches
// this only when that fits (small/medium models, the common Rfam case) and falls back to the
// global kernel for very long models.
//
// Layout of dynamic shared `smem`:
//   [0 .. Kp*(M+1))                    emission table, loaded cooperatively, shared by all threads
//   [Kp*(M+1) ..)                      DP rows, **interleaved** as dp[k*blockDim + tid] so a warp's
//                                      threads touch consecutive addresses for a fixed k (avoids
//                                      the bank conflicts a per-thread contiguous row would cause).
extern "C" __global__ void msv_u8_batch_kernel_smem(
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
    float*         __restrict__ out)
{
    extern __shared__ uint8_t smem[];
    uint8_t* s_rbv = smem;                          // Kp*(M+1)
    uint8_t* s_dp  = smem + (size_t)Kp * (M + 1);   // blockDim.x*(M+1), interleaved
    int      tid   = threadIdx.x;
    int      bd    = blockDim.x;
    int      j     = blockIdx.x * bd + tid;

    // Cooperatively stage the emission table (all threads, incl. out-of-range, participate).
    for (int idx = tid; idx < Kp * (M + 1); idx += bd) s_rbv[idx] = rbv[idx];
    __syncthreads();
    if (j >= n) return;   // no __syncthreads after this point, so early exit is safe

    int L    = lengths[j];
    int base = offsets[j];
    int tjb  = tjb_byte(scale_b, L);
    int tjbm = (tjb + tbm_b) & 0xff;

    for (int k = 0; k <= M; ++k) s_dp[k * bd + tid] = 0;

    int xJ = 0;
    int xB = ssub_u8(base_b, tjbm);

    for (int i = 1; i <= L; ++i) {
        int            x   = (int)dsq[base + (i - 1)];
        const uint8_t* row = s_rbv + (size_t)x * (M + 1);
        int            xE  = 0;

        for (int k = M; k >= 1; --k) {
            int prev = s_dp[(k - 1) * bd + tid];    // previous row's k-1 (not yet overwritten)
            int sv   = prev > xB ? prev : xB;
            sv = sadd_u8(sv, bias_b);
            sv = ssub_u8(sv, (int)row[k]);
            s_dp[k * bd + tid] = (uint8_t)sv;
            if (sv > xE) xE = sv;
        }

        if (sadd_u8(xE, bias_b) == 255) { out[j] = INFINITY; return; }

        int xJn = ssub_u8(xE, tec_b);
        if (xJn > xJ) xJ = xJn;
        int bmax = base_b > xJ ? base_b : xJ;
        xB = ssub_u8(bmax, tjbm);
    }

    out[j] = (((float)xJ - (float)tjb) - (float)base_b) / scale_b - 3.0f;
}

// ---- Viterbi filter (f32 max-plus, full M/I/D model) ----------------------------------------
// The middle cascade stage: the best single gapped local parse (`p7_GViterbi`), a faithful
// transcription of the scalar `ornament-hmm/src/fwd.rs::viterbi_nats` (max-plus semiring) — same
// transitions, same Dk→E folding, same length model — kept in f32 so it matches the scalar oracle
// to rounding. One thread per window (like MSV); it runs only on MSV survivors, so its aggregate
// volume is far lower.
//
// Memory: unlike MSV's single rolling row, Viterbi needs the previous *and* current rows of all
// three states — M(i,k) reads M/I/D(i-1,k-1), I(i,k) reads M/I(i-1,k), D(i,k) reads M/D(i,k-1) —
// so each thread keeps 6 rows of (M+1) floats in global scratch (stride 6*(M+1) per thread). This
// is the correct-first global-memory version; a shared-memory Viterbi (6× MSV's footprint) is a
// separate optimization.
//
// tsc: (M+1)*8 transition log-odds, tsc[k*8+type] with type MM=0,IM=1,DM=2,BM=3,MD=4,DD=5,MI=6,
// II=7 (HMMER p7P order). msc: [Kp][M+1] match-emission log-odds. Insert emissions are 0 (profile
// hardwires them). Length model (N/C/J loop/move) is recomputed on-device from L and nj; E move =
// E loop = ln(1/2) (multihit local). out[j] = raw Viterbi nats (host null-corrects to bits).

__device__ __forceinline__ float vmax(float a, float b) { return a >= b ? a : b; }

extern "C" __global__ void viterbi_batch_kernel(
    const float* __restrict__ tsc,   // (M+1)*8
    const float* __restrict__ msc,   // Kp*(M+1)
    int                        Kp,
    int                        M,
    float                      nj,
    const uint8_t* __restrict__ dsq,
    const int*     __restrict__ offsets,
    const int*     __restrict__ lengths,
    int                         n,
    float*         __restrict__ scratch,   // 6*(M+1) floats per thread
    float*         __restrict__ out)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;

    int   L    = lengths[j];
    int   base = offsets[j];
    int   S    = M + 1;
    float* buf = scratch + (size_t)j * 6 * S;
    float* Mp = buf;             // prev row M
    float* Ip = buf + S;         // prev row I
    float* Dp = buf + 2 * S;     // prev row D
    float* Mc = buf + 3 * S;     // cur row M
    float* Ic = buf + 4 * S;     // cur row I
    float* Dc = buf + 5 * S;     // cur row D

    // Length model (p7_ReconfigLength) + multihit-local E exits.
    float pmove = (2.0f + nj) / ((float)L + 2.0f + nj);
    float nloop = logf(1.0f - pmove);   // N/C/J loop
    float nmove = logf(pmove);          // N/C/J move
    float emove = logf(0.5f);           // E->C move  (= E loop)
    float eloop = emove;                // E->J loop

    // Row 0: M/I/D = -inf; specials.
    for (int k = 0; k <= M; ++k) { Mp[k] = NINF; Ip[k] = NINF; Dp[k] = NINF; }
    float xN = 0.0f, xJ = NINF, xC = NINF;
    float xB = nmove;   // S->N->B

    for (int i = 1; i <= L; ++i) {
        int          x  = (int)dsq[base + (i - 1)];
        const float* ms = msc + (size_t)x * S;
        float        xE = NINF;
        Mc[0] = NINF; Ic[0] = NINF; Dc[0] = NINF;

        for (int k = 1; k <= M; ++k) {
            const float* t  = tsc + (size_t)(k - 1) * 8; // MM,IM,DM,BM indexed by k-1
            const float* tk = tsc + (size_t)k * 8;       // MI,II indexed by k
            // match: max(M[k-1]+MM, I[k-1]+IM, B+BM, D[k-1]+DM) + emission
            float mv = vmax(vmax(Mp[k - 1] + t[0], Ip[k - 1] + t[1]),
                            vmax(xB + t[3], Dp[k - 1] + t[2]));
            Mc[k] = mv + ms[k];
            // insert (isc = 0), impossible at k = M
            Ic[k] = (k < M) ? vmax(Mp[k] + tk[6], Ip[k] + tk[7]) : NINF;
            // delete: max(Mcur[k-1]+MD, Dcur[k-1]+DD)  (same row, ascending k)
            Dc[k] = vmax(Mc[k - 1] + t[4], Dc[k - 1] + t[5]);
            // E gets Mk and Dk exits (esc = 0)
            xE = vmax(xE, vmax(Mc[k], Dc[k]));
        }

        // specials: J/C/N use prev-row values + this row's E; B uses this row's N and J.
        float xJn = vmax(xJ + nloop, xE + eloop);
        float xCn = vmax(xC + nloop, xE + emove);
        float xNn = xN + nloop;
        float xBn = vmax(xNn + nmove, xJn + nmove);
        xJ = xJn; xC = xCn; xN = xNn; xB = xBn;

        // swap prev/cur rows for the next iteration
        float* tmp;
        tmp = Mp; Mp = Mc; Mc = tmp;
        tmp = Ip; Ip = Ic; Ic = tmp;
        tmp = Dp; Dp = Dc; Dc = tmp;
    }

    out[j] = xC + nmove;   // + C->move; L == 0 ⇒ xC = -inf
}

// Shared-memory Viterbi. Same recurrence/result as `viterbi_batch_kernel`, but the per-thread DP
// state and the (block-shared) transition + emission tables live in shared memory. Profiling showed
// the global version is L2-latency-bound (2.8% compute, 42% L2, 167× slower than shared-mem MSV);
// this moves the hot dependent-load traffic into ~30-cycle shared memory.
//
// Two tricks keep the footprint down to **3** rows/thread (not the naive 6 = M/I/D × prev/cur):
//   * scan M,I **descending** in k so a single row per state doubles as prev (k-1, not yet
//     overwritten) and cur (k, written after reading the prev value at index k) — the MSV trick,
//     extended to I by reading M[k]/I[k] *before* overwriting them;
//   * compute D in a **second, ascending** pass (D(i,k) needs the current row's M[k-1]/D[k-1]),
//     after M/I are done — D stays at the previous row's values through pass 1 (where the M
//     recurrence reads Dprev[k-1]), then pass 2 overwrites it in place.
// So M/I/D are three single (M+1) rows updated in place; no per-row swap, and cur automatically
// becomes prev for the next row. Dynamic shared layout (floats):
//   [0 .. 8*S)              tsc          (S = M+1)
//   [8*S .. 8*S+Kp*S)       msc
//   [rest]                  dp, interleaved dp[(state*S + k)*blockDim + tid]  (conflict-free)
extern "C" __global__ void viterbi_batch_kernel_smem(
    const float* __restrict__ tsc,
    const float* __restrict__ msc_g,
    int                        Kp,
    int                        M,
    float                      nj,
    const uint8_t* __restrict__ dsq,
    const int*     __restrict__ offsets,
    const int*     __restrict__ lengths,
    int                         n,
    float*         __restrict__ out)
{
    extern __shared__ float sh[];
    int S   = M + 1;
    int bd  = blockDim.x;
    int tid = threadIdx.x;
    float* s_tsc = sh;
    float* s_msc = s_tsc + (size_t)S * 8;
    float* s_dp  = s_msc + (size_t)Kp * S;   // bd*3*S, interleaved

    for (int idx = tid; idx < S * 8; idx += bd) s_tsc[idx] = tsc[idx];
    for (int idx = tid; idx < Kp * S; idx += bd) s_msc[idx] = msc_g[idx];
    __syncthreads();

    int j = blockIdx.x * bd + tid;
    if (j >= n) return;

#define VM_(k) s_dp[((0 * S) + (k)) * bd + tid]
#define VI_(k) s_dp[((1 * S) + (k)) * bd + tid]
#define VD_(k) s_dp[((2 * S) + (k)) * bd + tid]

    int   L    = lengths[j];
    int   base = offsets[j];
    float pmove = (2.0f + nj) / ((float)L + 2.0f + nj);
    float nloop = logf(1.0f - pmove);
    float nmove = logf(pmove);
    float emove = logf(0.5f);
    float eloop = emove;

    for (int k = 0; k <= M; ++k) { VM_(k) = NINF; VI_(k) = NINF; VD_(k) = NINF; }  // row 0
    float xN = 0.0f, xJ = NINF, xC = NINF, xB = nmove;

    for (int i = 1; i <= L; ++i) {
        int          x  = (int)dsq[base + (i - 1)];
        const float* ms = s_msc + (size_t)x * S;
        float        xE = NINF;

        // Pass 1 (descending): cur M, I in place; D still holds the previous row (Dprev).
        for (int k = M; k >= 1; --k) {
            const float* t = s_tsc + (size_t)(k - 1) * 8;   // MM,IM,DM,BM by k-1
            float mv = vmax(vmax(VM_(k - 1) + t[0], VI_(k - 1) + t[1]),
                            vmax(xB + t[3], VD_(k - 1) + t[2]));
            mv += ms[k];
            float iv = NINF;
            if (k < M) {
                const float* tk = s_tsc + (size_t)k * 8;    // MI,II by k
                iv = vmax(VM_(k) + tk[6], VI_(k) + tk[7]);  // reads Mprev[k]/Iprev[k] before write
            }
            VM_(k) = mv;
            VI_(k) = iv;
            xE = vmax(xE, mv);
        }

        // Pass 2 (ascending): cur D in place. D(i,k) = max(Mcur[k-1]+MD, Dcur[k-1]+DD).
        VD_(0) = NINF;
        for (int k = 1; k <= M; ++k) {
            const float* t = s_tsc + (size_t)(k - 1) * 8;
            float dv = vmax(VM_(k - 1) + t[4], VD_(k - 1) + t[5]);
            VD_(k) = dv;
            xE = vmax(xE, dv);
        }

        float xJn = vmax(xJ + nloop, xE + eloop);
        float xCn = vmax(xC + nloop, xE + emove);
        float xNn = xN + nloop;
        float xBn = vmax(xNn + nmove, xJn + nmove);
        xJ = xJn; xC = xCn; xN = xNn; xB = xBn;
    }

    out[j] = xC + nmove;
#undef VM_
#undef VI_
#undef VD_
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

// Counter-free kernel diagnostics: registers/thread, static + dynamic shared bytes, and the
// theoretical occupancy (active blocks/SM × blockSize ÷ maxThreadsPerSM) from the CUDA occupancy
// API. This is a substitute for Nsight Compute's Occupancy section on clusters where GPU
// performance counters are locked (ERR_NVGPUCTRPERM). `Kp`/`M` size the u8 kernel's dynamic shared
// as `(Kp+blockDim)*(M+1)`. Prints one line per kernel to stdout.
static void ornament_gpu_kernel_line(const char* name, const void* fn, int block, size_t dyn_smem) {
    cudaFuncAttributes a;
    if (cudaFuncGetAttributes(&a, fn) != cudaSuccess) {
        printf("  %-26s <attr query failed>\n", name);
        cudaGetLastError();
        return;
    }
    int max_blocks = 0;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&max_blocks, fn, block, dyn_smem);
    int dev = 0;
    cudaGetDevice(&dev);
    int tpm = 0;
    cudaDeviceGetAttribute(&tpm, cudaDevAttrMaxThreadsPerMultiProcessor, dev);
    double occ = tpm > 0 ? (double)(max_blocks * block) / (double)tpm : 0.0;
    printf("  %-26s regs/thr=%-4d static_smem=%-6zu dyn_smem=%-7zu block=%-4d active_blk/SM=%-2d theo_occ=%.0f%%\n",
           name, a.numRegs, (size_t)a.sharedSizeBytes, dyn_smem, block, max_blocks, occ * 100.0);
    cudaGetLastError();
}

extern "C" void ornament_gpu_print_kernel_info(int Kp, int M) {
    size_t u8_smem = (size_t)(Kp + 128) * (M + 1);   // (Kp+block)*(M+1), block=128
    printf("kernel diagnostics (Kp=%d, M=%d):\n", Kp, M);
    ornament_gpu_kernel_line("msv_u8_batch_kernel_smem", (const void*)msv_u8_batch_kernel_smem, 128, u8_smem);
    ornament_gpu_kernel_line("msv_u8_batch_kernel",      (const void*)msv_u8_batch_kernel,      128, 0);
    ornament_gpu_kernel_line("msv_batch_kernel(f32)",    (const void*)msv_batch_kernel,         128, 0);
    ornament_gpu_kernel_line("viterbi_batch_kernel",     (const void*)viterbi_batch_kernel,     128, 0);
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

// ---- Streamed, double-buffered scoring over a resident strand -------------------------------
// The tile batch is split into fixed-size chunks and pipelined across S=2 CUDA streams so H2D
// (chunk i+1's offsets), kernel (chunk i), and D2H (chunk i-1's scores) overlap. Per-stream device
// + pinned host buffers are the double buffer; pinned staging is what lets cudaMemcpyAsync actually
// run concurrently with compute. The emission table (msc/rbv) is uploaded once by the caller and
// shared read-only by every chunk. Degrades cleanly to one chunk / one stream for small batches.
//
// `launch` enqueues the DP kernel for one chunk on a given stream — it captures the (already
// resident) emission table + strand + profile scalars, so this helper is agnostic to f32 vs uint8.

#define MAX_STREAMS      8
static const int DEFAULT_STREAMS = 2;
static const int DEFAULT_CHUNK   = 1 << 16;   // tiles/chunk (~64k); measured sweet spot on A30
                                              // (8k adds launch overhead; >=64k saturates)

// Read an int env var, clamped to [lo,hi]; `def` if unset/unparseable. Lets a benchmark A/B the
// stream count (ORNAMENT_GPU_STREAMS) and chunk size (ORNAMENT_GPU_CHUNK) without recompiling —
// STREAMS=1 is the no-overlap baseline, STREAMS>=2 is the double-buffered pipeline.
static int env_int(const char* name, int def, int lo, int hi) {
    const char* v = getenv(name);
    if (!v || !*v) return def;
    int x = atoi(v);
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return x;
}

// A reusable device context: streams + per-stream device/pinned buffers + a growable emission-
// table buffer, created ONCE and reused across many scoring calls (many models, many strands).
// A per-call scan (the old wrappers) re-created all of this every call — cheap once, but on an
// all-model scan (4k Rfam models) the stream/pinned setup dominated. The context hoists it out.
// off/len/out buffers are sized to `chunk` at creation; the scratch (global-kernel fallback) and
// emission buffers grow on demand as larger models appear.
struct GpuCtx {
    int          streams;
    int          chunk;
    cudaStream_t stream[MAX_STREAMS];
    int*         d_off[MAX_STREAMS];
    int*         d_len[MAX_STREAMS];
    float*       d_out[MAX_STREAMS];
    void*        d_scr[MAX_STREAMS];
    size_t       scr_cap[MAX_STREAMS];   // bytes
    int*         h_off[MAX_STREAMS];
    int*         h_len[MAX_STREAMS];
    float*       h_out[MAX_STREAMS];
    void*        d_emit;                 // reusable emission table (msc/rbv), grows on demand
    size_t       emit_cap;               // bytes
    void*        d_tsc;                  // reusable transition table (Viterbi), grows on demand
    size_t       tsc_cap;                // bytes
    int          smem_optin;             // opted-in per-block dynamic-shared max (bytes) for the u8 kernel
};

// Grow a device buffer to at least `need` bytes (free + realloc only when it must grow).
static int ensure_dev(void** p, size_t* cap, size_t need) {
    if (*cap >= need) return 0;
    if (*p) cudaFree(*p);
    *p = nullptr;
    *cap = 0;
    cudaError_t e = cudaMalloc(p, need);
    if (e != cudaSuccess) return (int)e;
    *cap = need;
    return 0;
}

static void ctx_free_internal(GpuCtx* c) {
    if (!c) return;
    for (int s = 0; s < MAX_STREAMS; ++s) {
        if (c->d_off[s]) cudaFree(c->d_off[s]);
        if (c->d_len[s]) cudaFree(c->d_len[s]);
        if (c->d_out[s]) cudaFree(c->d_out[s]);
        if (c->d_scr[s]) cudaFree(c->d_scr[s]);
        if (c->h_off[s]) cudaFreeHost(c->h_off[s]);
        if (c->h_len[s]) cudaFreeHost(c->h_len[s]);
        if (c->h_out[s]) cudaFreeHost(c->h_out[s]);
        if (c->stream[s]) cudaStreamDestroy(c->stream[s]);
    }
    if (c->d_emit) cudaFree(c->d_emit);
    if (c->d_tsc) cudaFree(c->d_tsc);
    free(c);
}

// Create a context. streams/chunk <= 0 fall back to ORNAMENT_GPU_STREAMS/CHUNK (else defaults).
extern "C" int ornament_gpu_ctx_create(int streams, int chunk, void** out) {
    *out = nullptr;
    if (streams <= 0) streams = env_int("ORNAMENT_GPU_STREAMS", DEFAULT_STREAMS, 1, MAX_STREAMS);
    if (streams > MAX_STREAMS) streams = MAX_STREAMS;
    if (chunk <= 0) chunk = env_int("ORNAMENT_GPU_CHUNK", DEFAULT_CHUNK, 1, 1 << 24);

    GpuCtx* c = (GpuCtx*)calloc(1, sizeof(GpuCtx));
    if (!c) return (int)cudaErrorMemoryAllocation;
    c->streams = streams;
    c->chunk = chunk;

    int rc = 0;
    for (int s = 0; s < streams; ++s) {
        if ((rc = (int)cudaStreamCreate(&c->stream[s]))) goto fail;
        if ((rc = (int)cudaMalloc(&c->d_off[s], (size_t)chunk * sizeof(int)))) goto fail;
        if ((rc = (int)cudaMalloc(&c->d_len[s], (size_t)chunk * sizeof(int)))) goto fail;
        if ((rc = (int)cudaMalloc(&c->d_out[s], (size_t)chunk * sizeof(float)))) goto fail;
        if ((rc = (int)cudaHostAlloc(&c->h_off[s], (size_t)chunk * sizeof(int), cudaHostAllocDefault))) goto fail;
        if ((rc = (int)cudaHostAlloc(&c->h_len[s], (size_t)chunk * sizeof(int), cudaHostAllocDefault))) goto fail;
        if ((rc = (int)cudaHostAlloc(&c->h_out[s], (size_t)chunk * sizeof(float), cudaHostAllocDefault))) goto fail;
    }

    // Opt the shared-memory u8 kernel into the device's larger per-block dynamic-shared limit
    // (sm_80: ~164 KB vs the 48 KB static default), so longer models can still run on shared
    // memory. Best-effort: on any failure keep the safe 48 KB budget. Device 0 (single-GPU runs).
    {
        int optin = 48 * 1024;
        int dev = 0;
        cudaGetDevice(&dev);
        int attr = 0;
        if (cudaDeviceGetAttribute(&attr, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev)
                == cudaSuccess
            && attr > optin
            && cudaFuncSetAttribute((const void*)msv_u8_batch_kernel_smem,
                                    cudaFuncAttributeMaxDynamicSharedMemorySize, attr)
                   == cudaSuccess) {
            optin = attr;
            // Same opt-in for the shared-memory Viterbi kernel (best-effort).
            cudaFuncSetAttribute((const void*)viterbi_batch_kernel_smem,
                                 cudaFuncAttributeMaxDynamicSharedMemorySize, attr);
        } else {
            cudaGetLastError(); // clear any sticky error from the probe
        }
        c->smem_optin = optin;
    }

    *out = c;
    return 0;
fail:
    ctx_free_internal(c);
    return rc;
}

extern "C" int ornament_gpu_ctx_destroy(void* ctx) {
    ctx_free_internal((GpuCtx*)ctx);
    return 0;
}

// Streamed, double-buffered scoring using a context's buffers (no per-call alloc/free). The tile
// batch is chunked and round-robined across the context's streams so chunk i+1's H2D, chunk i's
// kernel, and chunk i-1's D2H overlap. `launch` enqueues the DP kernel for one chunk (agnostic to
// f32 vs uint8; it captures the resident strand + emission table + profile scalars).
template <typename Launch>
static int ctx_stream_run(
    GpuCtx*      ctx,
    int          n,
    int          M,
    size_t       elem_size,   // scratch bytes per DP cell (global-kernel fallback only)
    int          threads,     // block size (chosen by the caller: 128 for f32/global, ≤128 for u8 smem)
    const int*   starts,
    const int*   lengths,
    float*       out,
    Launch       launch)
{
    int streams = ctx->streams;
    int chunk = (n < ctx->chunk) ? n : ctx->chunk;

    size_t scratch_bytes = (size_t)chunk * (M + 1) * elem_size;
    for (int s = 0; s < streams; ++s) {
        int rc = ensure_dev(&ctx->d_scr[s], &ctx->scr_cap[s], scratch_bytes);
        if (rc) return rc;
    }

    int  pend_off[MAX_STREAMS];
    int  pend_n[MAX_STREAMS];
    bool pend[MAX_STREAMS] = {false};

    for (int off = 0; off < n; off += chunk) {
        int s = (off / chunk) % streams;
        if (pend[s]) {
            cudaError_t e = cudaStreamSynchronize(ctx->stream[s]);
            if (e != cudaSuccess) return (int)e;
            memcpy(out + pend_off[s], ctx->h_out[s], (size_t)pend_n[s] * sizeof(float));
            pend[s] = false;
        }
        int cn = (n - off < chunk) ? (n - off) : chunk;
        memcpy(ctx->h_off[s], starts + off,  (size_t)cn * sizeof(int));
        memcpy(ctx->h_len[s], lengths + off, (size_t)cn * sizeof(int));
        cudaError_t e;
        if ((e = cudaMemcpyAsync(ctx->d_off[s], ctx->h_off[s], (size_t)cn * sizeof(int),
                                 cudaMemcpyHostToDevice, ctx->stream[s]))) return (int)e;
        if ((e = cudaMemcpyAsync(ctx->d_len[s], ctx->h_len[s], (size_t)cn * sizeof(int),
                                 cudaMemcpyHostToDevice, ctx->stream[s]))) return (int)e;
        int blocks = (cn + threads - 1) / threads;
        launch(ctx->stream[s], blocks, threads, ctx->d_off[s], ctx->d_len[s], cn,
               ctx->d_scr[s], ctx->d_out[s]);
        if ((e = cudaGetLastError())) return (int)e;
        if ((e = cudaMemcpyAsync(ctx->h_out[s], ctx->d_out[s], (size_t)cn * sizeof(float),
                                 cudaMemcpyDeviceToHost, ctx->stream[s]))) return (int)e;
        pend_off[s] = off;
        pend_n[s]   = cn;
        pend[s]     = true;
    }

    for (int s = 0; s < streams; ++s) {
        if (pend[s]) {
            cudaError_t e = cudaStreamSynchronize(ctx->stream[s]);
            if (e != cudaSuccess) return (int)e;
            memcpy(out + pend_off[s], ctx->h_out[s], (size_t)pend_n[s] * sizeof(float));
        }
    }
    return 0;
}

// f32 MSV over a resident strand, using a reusable context. msc: Kp*(M+1) floats (host); uploaded
// into the context's growable emission buffer. `starts` are 0-based strand offsets, `lengths` the
// window lengths.
extern "C" int ornament_gpu_ctx_msv_f32(
    void*          ctx_,
    const float*   msc,
    int            Kp,
    int            M,
    const void*    d_strand,
    const int*     starts,
    const int*     lengths,
    int            n,
    float*         out)
{
    if (n <= 0) return 0;
    GpuCtx* ctx = (GpuCtx*)ctx_;
    size_t msc_bytes = (size_t)Kp * (M + 1) * sizeof(float);
    int rc = ensure_dev(&ctx->d_emit, &ctx->emit_cap, msc_bytes);
    if (rc) return rc;
    cudaError_t e = cudaMemcpy(ctx->d_emit, msc, msc_bytes, cudaMemcpyHostToDevice);
    if (e != cudaSuccess) return (int)e;
    const float* d_msc = (const float*)ctx->d_emit;

    return ctx_stream_run(
        ctx, n, M, sizeof(float), 128, starts, lengths, out,
        [=](cudaStream_t st, int blocks, int threads,
            const int* d_off, const int* d_len, int cn, void* d_scr, float* d_out) {
            msv_batch_kernel<<<blocks, threads, 0, st>>>(
                d_msc, Kp, M, (const uint8_t*)d_strand, d_off, d_len, cn,
                (float*)d_scr, d_out);
        });
}

// Pick the block size for the shared-memory uint8 kernel: the largest of {128,64,32} whose
// footprint `(Kp+bd)*(M+1)` fits the opted-in shared limit. Smaller blocks cover longer models
// (fewer per-block DP rows) at the cost of occupancy — still far faster than the global kernel.
// Returns 0 if even 32 threads won't fit (caller falls back to the global kernel).
static int pick_u8_threads(int Kp, int M, size_t smem_limit) {
    for (int bd = 128; bd >= 32; bd >>= 1) {
        if ((size_t)(Kp + bd) * (M + 1) <= smem_limit) return bd;
    }
    return 0;
}

// Block size for the shared-memory Viterbi kernel: largest of {128,64,32} whose dynamic-shared
// footprint `(8 + Kp + 3*bd)*(M+1)` floats fits the opted-in limit (tsc + msc + 3 DP rows/thread).
// Viterbi's f32 3-row state is 12× the u8 MSV's 1-byte row, so it fits fewer/shorter models than
// MSV; 0 ⇒ fall back to the global Viterbi kernel.
static int pick_vit_threads(int Kp, int M, size_t smem_limit) {
    for (int bd = 128; bd >= 32; bd >>= 1) {
        size_t bytes = (size_t)(8 + Kp + 3 * bd) * (M + 1) * sizeof(float);
        if (bytes <= smem_limit) return bd;
    }
    return 0;
}

// Reduced-precision (uint8) MSV over a resident strand, using a reusable context. rbv: Kp*(M+1)
// bytes (host). Prefers the shared-memory kernel when it fits (ORNAMENT_GPU_SMEM=0 forces global).
extern "C" int ornament_gpu_ctx_msv_u8(
    void*          ctx_,
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
    if (n <= 0) return 0;
    GpuCtx* ctx = (GpuCtx*)ctx_;
    size_t rbv_bytes = (size_t)Kp * (M + 1);
    int rc = ensure_dev(&ctx->d_emit, &ctx->emit_cap, rbv_bytes);
    if (rc) return rc;
    cudaError_t e = cudaMemcpy(ctx->d_emit, rbv, rbv_bytes, cudaMemcpyHostToDevice);
    if (e != cudaSuccess) return (int)e;
    const uint8_t* d_rbv = (const uint8_t*)ctx->d_emit;

    int use_smem = env_int("ORNAMENT_GPU_SMEM", 1, 0, 1);
    // Shared-mem budget = the context's opted-in per-block max (sm_80: ~164 KB vs 48 KB default),
    // set once at ctx_create. Choose the largest block size that fits; 0 ⇒ use the global kernel.
    size_t smem_limit = (size_t)ctx->smem_optin;
    int smem_threads = use_smem ? pick_u8_threads(Kp, M, smem_limit) : 0;
    int threads = smem_threads ? smem_threads : 128;

    return ctx_stream_run(
        ctx, n, M, sizeof(uint8_t), threads, starts, lengths, out,
        [=](cudaStream_t st, int blocks, int th,
            const int* d_off, const int* d_len, int cn, void* d_scr, float* d_out) {
            if (smem_threads) {
                size_t shmem = (size_t)(Kp + th) * (M + 1);
                msv_u8_batch_kernel_smem<<<blocks, th, shmem, st>>>(
                    d_rbv, Kp, M, scale_b, base_b, bias_b, tbm_b, tec_b,
                    (const uint8_t*)d_strand, d_off, d_len, cn, d_out);
            } else {
                msv_u8_batch_kernel<<<blocks, th, 0, st>>>(
                    d_rbv, Kp, M, scale_b, base_b, bias_b, tbm_b, tec_b,
                    (const uint8_t*)d_strand, d_off, d_len, cn, (uint8_t*)d_scr, d_out);
            }
        });
}

// Per-call wrappers (no persistent context): spin up a temporary context, score, tear it down.
// Kept for the simple one-shot API and the parity tests; a many-model scan should reuse a context.
extern "C" int ornament_gpu_msv_batch_resident(
    const float* msc, int Kp, int M, const void* d_strand,
    const int* starts, const int* lengths, int n, float* out)
{
    if (n <= 0) return 0;
    void* ctx = nullptr;
    int rc = ornament_gpu_ctx_create(0, 0, &ctx);
    if (rc) return rc;
    rc = ornament_gpu_ctx_msv_f32(ctx, msc, Kp, M, d_strand, starts, lengths, n, out);
    ornament_gpu_ctx_destroy(ctx);
    return rc;
}

extern "C" int ornament_gpu_msv_u8_batch_resident(
    const uint8_t* rbv, int Kp, int M, float scale_b, int base_b, int bias_b, int tbm_b, int tec_b,
    const void* d_strand, const int* starts, const int* lengths, int n, float* out)
{
    if (n <= 0) return 0;
    void* ctx = nullptr;
    int rc = ornament_gpu_ctx_create(0, 0, &ctx);
    if (rc) return rc;
    rc = ornament_gpu_ctx_msv_u8(ctx, rbv, Kp, M, scale_b, base_b, bias_b, tbm_b, tec_b,
                                 d_strand, starts, lengths, n, out);
    ornament_gpu_ctx_destroy(ctx);
    return rc;
}

// Viterbi (f32) over a resident strand, using a reusable context. tsc: (M+1)*8 transition
// log-odds; msc: Kp*(M+1) match-emission log-odds; nj: expected-hit count (length model). Both
// tables are uploaded into growable context buffers; scratch is 6*(M+1) floats per tile.
extern "C" int ornament_gpu_ctx_vit_f32(
    void*          ctx_,
    const float*   tsc,
    const float*   msc,
    int            Kp,
    int            M,
    float          nj,
    const void*    d_strand,
    const int*     starts,
    const int*     lengths,
    int            n,
    float*         out)
{
    if (n <= 0) return 0;
    GpuCtx* ctx = (GpuCtx*)ctx_;
    size_t tsc_bytes = (size_t)(M + 1) * 8 * sizeof(float);
    size_t msc_bytes = (size_t)Kp * (M + 1) * sizeof(float);
    int rc = ensure_dev(&ctx->d_tsc, &ctx->tsc_cap, tsc_bytes);
    if (rc) return rc;
    rc = ensure_dev(&ctx->d_emit, &ctx->emit_cap, msc_bytes);
    if (rc) return rc;
    cudaError_t e;
    if ((e = cudaMemcpy(ctx->d_tsc, tsc, tsc_bytes, cudaMemcpyHostToDevice))) return (int)e;
    if ((e = cudaMemcpy(ctx->d_emit, msc, msc_bytes, cudaMemcpyHostToDevice))) return (int)e;
    const float* d_tsc = (const float*)ctx->d_tsc;
    const float* d_msc = (const float*)ctx->d_emit;

    // Prefer the shared-memory kernel when its footprint fits; else the global-scratch kernel.
    // ORNAMENT_GPU_SMEM=0 forces global (for A/B). Global path needs 6*(M+1) floats of scratch per
    // tile; the shared path needs none (elem_size 0 ⇒ ctx_stream_run skips scratch alloc).
    int use_smem = env_int("ORNAMENT_GPU_SMEM", 1, 0, 1);
    int vit_threads = use_smem ? pick_vit_threads(Kp, M, (size_t)ctx->smem_optin) : 0;
    int threads = vit_threads ? vit_threads : 128;
    size_t elem = vit_threads ? 0 : 6 * sizeof(float);

    return ctx_stream_run(
        ctx, n, M, elem, threads, starts, lengths, out,
        [=](cudaStream_t st, int blocks, int th,
            const int* d_off, const int* d_len, int cn, void* d_scr, float* d_out) {
            if (vit_threads) {
                size_t shmem = (size_t)(8 + Kp + 3 * th) * (M + 1) * sizeof(float);
                viterbi_batch_kernel_smem<<<blocks, th, shmem, st>>>(
                    d_tsc, d_msc, Kp, M, nj, (const uint8_t*)d_strand, d_off, d_len, cn, d_out);
            } else {
                viterbi_batch_kernel<<<blocks, th, 0, st>>>(
                    d_tsc, d_msc, Kp, M, nj, (const uint8_t*)d_strand, d_off, d_len, cn,
                    (float*)d_scr, d_out);
            }
        });
}

extern "C" int ornament_gpu_vit_batch_resident(
    const float* tsc, const float* msc, int Kp, int M, float nj,
    const void* d_strand, const int* starts, const int* lengths, int n, float* out)
{
    if (n <= 0) return 0;
    void* ctx = nullptr;
    int rc = ornament_gpu_ctx_create(0, 0, &ctx);
    if (rc) return rc;
    rc = ornament_gpu_ctx_vit_f32(ctx, tsc, msc, Kp, M, nj, d_strand, starts, lengths, n, out);
    ornament_gpu_ctx_destroy(ctx);
    return rc;
}
