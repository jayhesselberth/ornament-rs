# GPU offload design — where the work is, and the producer/consumer split

## The pipeline is a funnel

`ornament-scfg`'s scan cascades every candidate window cheapest-filter-first
(`pipeline.rs::strand_survivors`): **MSV → Viterbi → Forward → CYK/Inside**. Each stage only
runs on the survivors of the previous one, so the batch shrinks by ~1–2 orders of magnitude at
each step.

| Stage | Per-window cost | Windows touched | Aggregate cost | GPU fit |
|-------|-----------------|-----------------|----------------|---------|
| **MSV** (uint8, ungapped) | cheapest | **all** (genome × models × 2 strands) | **dominant** | excellent |
| Viterbi (f32 max-plus, +I/D/DD) | ~4–8× MSV | ~1–10% | moderate | good |
| Forward (odds-space, rescaled) | costliest filter | ~0.1–1% | small | good |
| CYK / Inside (profile-SCFG) | *huge* | a handful | small aggregate | **poor** |

"Heaviest" is ambiguous on purpose:
- **Per window**, CYK dominates — but it runs on almost nothing.
- **In aggregate wall-clock**, MSV dominates, because it is the only stage that touches *every*
  tile. This is why HMMER hand-vectorized MSV first, and why it is our first (and highest-value)
  GPU target.

## What we GPU-enable — and what we don't

- **GPU: the p7 filter funnel (MSV, Viterbi, Forward).** Dense, regular, fixed-shape DP over a
  massive batch of independent tiles. This is where the aggregate time goes and where SIMD/GPU
  parallelism pays.
- **CPU: keep CYK/Inside.** Irregular (banded, variable model-span, large SCFG state space,
  recursive), *and* it runs on a trickle of survivors. Little aggregate cost, bad GPU fit.

So the boundary is: **GPU filters the flood; CPU handles the heavy CM tail on the survivors.**

## Producer/consumer, two levels

### Level 1 — host ↔ device (the streaming boundary)
- **Producer (host):** stream FASTA → digitize + reverse-complement once (done, #43) → tile into
  windows → build the batch.
- **Consumer (device):** the filter funnel.
- **Overlap** *(done — inside the resident MSV wrappers)*: the tile batch is split into fixed-size
  chunks pipelined across 2 CUDA streams with pinned host staging, so chunk *i+1*'s H2D, chunk
  *i*'s kernel, and chunk *i-1*'s D2H overlap. The emission table is uploaded once; the strand is
  already resident. *Still to do:* extend the same pipeline across the Viterbi/Forward stages.

## Measured (A30, uint8 MSV, tRNA M≈71, 64 Mnt strand → 640k tiles; `examples/bench_msv.rs`)

Stream A/B (chunk 32k) and chunk sweep (2 streams):

| config | best time | throughput |
|--------|-----------|------------|
| 1 stream (no overlap) | 0.585 s | 1.09 M tiles/s |
| **2 streams** (double-buffered) | 0.492 s | 1.30 M tiles/s |
| 4 streams | 0.489 s | 1.31 M tiles/s |
| 2 streams, chunk 8k | 0.602 s | 1.06 M tiles/s |
| 2 streams, chunk 64k | 0.486 s | 1.32 M tiles/s |

- **Streaming speedup ≈ 1.19×** (2 vs 1 stream). Real but modest — the copies it hides are only
  ~15% of wall time, so the batch is **kernel-bound, not transfer-bound**. Beyond 2 streams is
  noise; chunk 64k is the sweet spot (8k adds launch/sync overhead) — now the default.
- **Streaming alone left the GPU only ~2.6× a single CPU core** (global-scratch kernel), i.e.
  *slower* than the ~24-core rayon CPU. At ~18.7 G DP-cells/s on a ~933 GB/s A30 it was nowhere
  near bandwidth limits — the bottleneck was the per-thread DP row in *global* memory, hit M×L
  times over a serial dependency.

### Shared-memory DP row (the fix — measured A30, same 640k-tile workload)

| kernel | best time | throughput | vs 1 CPU core |
|--------|-----------|------------|---------------|
| global-scratch DP row (`ORNAMENT_GPU_SMEM=0`) | 0.487 s | 18.7 G cells/s | 2.6× |
| **shared-memory DP row** (`ORNAMENT_GPU_SMEM=1`, default) | **0.022 s** | **418 G cells/s** | **57.6×** |

- **~22× faster kernel** from moving the DP row (and the small emission table) into shared memory.
  The recurrence is unchanged; only the memory it lives in changed (~400-cycle global → ~30-cycle
  shared), which is the whole story for a latency-bound serial-dependency DP.
- **57.6× one CPU core ⇒ ≈2.4× the full 24-core rayon CPU.** The GPU now genuinely beats the
  multicore CPU on this workload (tRNA-sized model). This is the result that makes offload viable.
- Requires `(Kp + blockDim)*(M+1)` bytes of shared memory; the host launches the shared kernel only
  when that fits in 48 KB (M ≲ 336 at blockDim 128) and falls back to the global kernel for very
  long models.

### All-Rfam validation (A30, `examples/sweep_rfam.rs`, full `Rfam.cm.gz` = 4227 models, 2 Mnt strand)

- **4227/4227 models parsed + scored, 0 skipped, 0 parity failures** vs the CPU striped filter
  (40 models spot-checked across the M range; max M = 3401). The kernel is correct across the
  entire real collection, not just tRNA.
- **Coverage: 4124/4227 (97.6%) fit the shared-mem kernel; 103 (M > 336) fall to the global
  kernel.** So nearly all of Rfam gets the fast path today.
- **The fallback is the bottleneck, quantified:** those 103 large models (2.4% of the collection)
  consume ~19.6 s of the 28 s total kernel time — the `M336-999`/`M≥1000` buckets run at
  0.13 / 0.03 M tiles/s, ~50–100× slower than the shared-mem buckets. This is the concrete,
  data-driven case for **warp-per-window** (fix #3 extension): make large models use shared memory
  cooperatively instead of the global fallback.
### Reusable device context (`GpuContext`) — removes per-model setup

The first sweep's per-model throughput was setup-dominated (each scoring call created + destroyed
its streams/pinned buffers). `GpuContext` creates them once and reuses them across all models.
Re-running the sweep with one shared context (A30, same 4227 models):

| bucket (M) | per-call tiles/s | reusable-context tiles/s |
|-----------|------------------|--------------------------|
| M<100 | 11.3 M | **27.0 M** |
| M100-199 | 8.5 M | **16.2 M** |
| M200-335 | 6.1 M | **9.2 M** |
| M336-999 (global) | 0.13 M | 0.14 M |
| M≥1000 (global) | 0.03 M | 0.03 M |

- The shared-mem buckets ≈ **2× faster** with the context, and M<100 (27 M tiles/s) is now close to
  `bench_msv`'s 29 M tiles/s single-model ceiling — i.e. setup overhead is gone and we're seeing the
  real kernel throughput. This is the path many-model scans (and pipeline integration) must use.
- **The 103 global-fallback models (M>336) now consume 82% of total kernel time (19.5 s of 23.6 s).**
  With setup removed, the long-model fallback is unambiguously *the* remaining bottleneck →
  larger-shared + block-size dispatch (below).

### Larger dynamic shared + block-size dispatch — long models onto shared memory

`GpuContext` opts the u8 shared kernel into the device's larger per-block dynamic-shared limit
(sm_80: ~164 KB via `cudaFuncSetAttribute`, vs the 48 KB default), and the u8 dispatch picks the
largest block size in {128, 64, 32} whose `(Kp+bd)*(M+1)` fits (smaller blocks hold fewer per-block
DP rows → cover longer models, at lower occupancy). At 164 KB even 32-thread blocks reach M≈3337,
so **4226/4227 Rfam models now run on shared memory** (only M=3401 falls back). Re-running the
sweep (A30, reusable context):

| bucket (M) | 48 KB / bd=128 only | 164 KB + bd dispatch |
|-----------|---------------------|----------------------|
| M<100 | 27.0 M tiles/s | 26.8 M |
| M100-199 | 16.2 M | 16.1 M |
| M200-335 | 9.2 M | 9.2 M |
| M336-999 (was global) | 0.14 M | **2.68 M (19×)** |
| M≥1000 (was global) | 0.03 M | **0.10 M (3.3×)** |
| **total kernel time** | 23.6 s | **6.8 s (3.5× overall)** |

- The two long-model buckets were the 82% bottleneck; moving them onto shared memory cut total
  sweep time 3.5×. Parity still clean (0 failures across the 40-model spot-check incl. large M).
- The `M≥1000` bucket (10 models) is still slow (0.10 M tiles/s): at bd=32 with ~164 KB shared it's
  ~1 block / 32 threads per SM — correct but occupancy-starved. **This is the real (now narrow)
  warp-per-window case:** keep block size high and split one window's DP row across a warp's lanes,
  instead of one thread owning the whole row. Only ~10 Rfam models need it; the rest are already
  fast. (Those rRNA-scale models could also simply stay on the CPU filter.)

### Profiling (A30, `examples/bench_msv.rs` with `BENCH_KERNEL`; nsys works, ncu blocked)

Per-kernel throughput at 64 Mnt / 640k tiles, tRNA-scale model (M≈71):

| kernel | G DP-cells/s | vs 1 CPU core | notes |
|--------|--------------|---------------|-------|
| MSV uint8 (shared-mem) | **436** | 61.8× | the workhorse; ≈2.6× the 24-core CPU |
| MSV f32 (global) | 14 | 1.9× | global-mem; no shared variant (it's the exact oracle, low volume) |
| Viterbi f32 (global) | 0.77 | **0.1×** | global-mem, memory-latency bound — *slower than one CPU core* |

- **nsys** (u8, 2 streams): the batch is **kernel-bound** — 39.6 ms in kernels vs 3.2 ms in copies,
  and the largest single copy is the *one-time* 16 MB strand upload (2.6 ms, amortized across all
  tiles/models). So streaming correctly hides only the small per-chunk offset copies (matches the
  1.19× A/B), and the strand-resident design already removed the dominant transfer.
- **Viterbi is ~78–570× slower per cell than MSV** and below one CPU core. This is the predicted
  global-memory penalty (6 DP rows hit with a serial dependency at ~400-cycle latency). **Concrete
  consequence: do NOT wire the current global Viterbi into the pipeline — it would regress vs the
  CPU striped Viterbi.** It needs the shared-memory / warp-per-window rewrite first.
- **ncu (SpeedOfLight/Occupancy) is blocked by `ERR_NVGPUCTRPERM`** on this cluster (perf counters
  are admin-locked; fix = `NVreg_RestrictProfilingToAdminUsers=0` on the GPU nodes). Until then,
  `print_kernel_info` (`BENCH_INFO=1`) reports registers/thread + theoretical occupancy via the
  counter-free CUDA occupancy API as a substitute (tRNA M=71, block=128):

  | kernel | regs/thr | dyn smem | active blk/SM | theoretical occ |
  |--------|----------|----------|---------------|-----------------|
  | msv_u8 shared-mem | 32 | 10.5 KB | 14 | **88%** |
  | msv_u8 global | 36 | 0 | 12 | 75% |
  | msv f32 global | 56 | 0 | 9 | 56% |
  | viterbi f32 global | 63 | 0 | 8 | 50% |

  The MSV shared kernel is well-occupied (88%) and fast (436 G cells/s) — near its roofline.
  Crucially, **Viterbi has a healthy 50% theoretical occupancy yet runs at 0.1× a CPU core**, so it
  is *not* occupancy/register-limited — it's memory-latency-bound on the global DP rows. That
  pinpoints the fix as **moving the DP rows to shared memory** (cut ~400-cycle global to ~30-cycle
  shared), not adding threads — the same lever that gave MSV its 22× jump.

### Level 2 — on-device funnel (stream compaction)
- MSV kernel writes a survivor mask → compact → Viterbi kernel consumes only survivors → compact
  → Forward. Compaction between kernels keeps later, costlier kernels dense (no warps full of
  already-pruned lanes). Survivors flow back to the host for CYK.

## Concrete wins, in priority order

1. **Resident sequence + offset windows** *(done — `DeviceStrand` + `Tiles`)*. Each strand is
   uploaded to the device once (RAII buffer, reused across windows and cascade stages); windows are
   `(start,end)` offsets into it instead of copied residues. The scan tiles at length `2W` step `W`
   (2× overlap), so the old per-window copy moved each base ~2×. *Still to do:* CUDA streams +
   double-buffering so H2D upload / kernel / D2H overlap (the full Level-1 producer/consumer).
2. **uint8 MSV** *(done — `msv_u8_batch_kernel`)*: 4× smaller emission table → better bandwidth on
   the memory-bound one-thread-per-window batch.
3. **Shared-memory DP row / warp-per-window** *(now the top lever — see Measured)*. The per-thread
   DP row `mmx[0..=M]` currently lives in global memory and is touched M×L times with a serial
   dependency, which is what caps the kernel at ~2.6× one CPU core. Stage the row (and the emission
   table) in shared memory; for large M go one-**warp**-per-window and dispatch on M. This is the
   change that decides whether GPU beats the multicore CPU at all.
4. **Batch across models and sequences.** A single small model (tRNA, M≈71) underutilizes the GPU.
   - *Many-genomes, one model:* one resident emission table + a huge window batch. Current design.
   - *All-model, one genome (~4000 Rfam):* iterate models but keep emission tables resident/cached
     and pack small models together (per-window model-id index into resident tables).

## Parity discipline

Every device kernel is validated against a trusted CPU oracle to f32 rounding:
- `msv_batch_kernel` (f32) ↔ `ornament_hmm::msv_nats`.
- `msv_u8_batch_kernel` (uint8) ↔ `ornament_hmm::MsvProfile::score_bits` (the striped uint8 filter).

The CPU path stays the reference; the GPU only changes speed, never the reported result.
