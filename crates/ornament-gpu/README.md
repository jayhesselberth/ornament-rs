# ornament-gpu

CUDA-accelerated p7 filter kernels — a **prototype** toward offloading the all-model /
many-genomes `scan` to GPU.

The scan's heaviest, most parallel workload is the p7 prefilter cascade run over a huge batch
of `(model, window)` tiles. This crate begins with the cheapest, most GPU-friendly stage — the
**MSV filter** — batched **one thread per window** on the device (`kernels/msv.cu`). The kernel
is a faithful f32 max-plus transcription of the scalar [`ornament_hmm::msv_nats`] oracle, so the
device result matches the CPU to f32 rounding (verified by the `gpu_matches_cpu` parity test).

See [`DESIGN.md`](DESIGN.md) for where the work is heaviest and the producer/consumer split.

## Status

- [x] Batch MSV kernel (f32, one-thread-per-window) + CUDA-Runtime-API host wrappers.
- [x] Feature-gated build (`nvcc` via `build.rs`), CPU oracle, GPU↔CPU parity test.
- [x] uint8 reduced-precision MSV (throughput; matches the striped `ornament-hmm::MsvProfile`).
- [x] Resident sequence + offset windows (`DeviceStrand` + `Tiles`; kills the 2× overlap copy).
- [x] CUDA streams + double-buffering: chunked tile batch pipelined across 2 streams with pinned
      staging so H2D / kernel / D2H overlap (internal to the resident wrappers).
- [x] Benchmark (`examples/bench_msv.rs`) + measured the streaming speedup (~1.19×) — see DESIGN.md.
      Finding: the batch is kernel-bound, and the kernel is only ~2.6× one CPU core.
- [x] **Shared-memory DP row** (`msv_u8_batch_kernel_smem`): ~22× faster than the global-scratch
      kernel (18.7 → 418 G DP-cells/s on an A30), ≈2.4× the full 24-core CPU. Default when it fits
      in 48 KB shared (M ≲ 336 @ blockDim 128); falls back to the global kernel for very long models.
      `ORNAMENT_GPU_SMEM=0` forces global for A/B.
- [x] **Validated across the full Rfam CM collection** (`examples/sweep_rfam.rs`): 4227/4227 models
      scored, 0 parity failures (max M=3401); 4124 (97.6%) fit shared-mem, 103 (M>336) fall to
      global — and those 103 eat ~70% of kernel time. See DESIGN.md.
- [x] Reusable device context (`GpuContext`): streams + buffers created once, reused across models
      — ~2× the shared-mem buckets in the Rfam sweep (M<100: 11→27 M tiles/s, ≈ single-model
      ceiling). Old per-call `*_batch_resident` are now thin temp-context shims. See DESIGN.md.
- [x] **Larger dynamic shared (~164 KB opt-in) + block-size dispatch {128/64/32}**: 4226/4227 Rfam
      models now run on shared memory (only M=3401 falls back). M336-999 bucket 19× faster, total
      Rfam sweep 3.5× faster (23.6→6.8 s). See DESIGN.md.
- [x] Wired `GpuContext` into the `ornament-scfg` scan pipeline as an optional backend (the
      `cuda` feature + `PipelineParams.gpu_msv` / `ORNAMENT_GPU=1`); GPU MSV pass-mask, CPU
      Vit/Fwd/CYK on survivors, hits identical. `gpu_msv_pipeline_matches_cpu` passes on A30.
- [x] Viterbi kernel (`viterbi_batch_kernel`): f32 max-plus M/I/D DP, one thread/window, global
      scratch (6 rows/thread); parity vs scalar `viterbi_nats` (`gpu_viterbi_matches_scalar`).
- [x] **Shared-memory Viterbi** (`viterbi_batch_kernel_smem`): 3 in-place DP rows (2-pass:
      descending M/I, ascending D) + tables in shared; blockDim dispatch {128/64/32}. **57× the
      global kernel** (0.86→49.4 G cells/s), 0.1×→**7.0× one CPU core**; same parity test. See DESIGN.md.
- [x] Tile aggregation (`DeviceStrand::upload_concat` + `Tiles::push_segment`): pack many
      sequences/strands into one resident buffer + one launch (no kernel change — offsets already
      work). Parity-tested (`gpu_aggregation_matches_separate`); `bench_agg` shows 178–469× over
      per-sequence launches for many small sequences (eliminates per-launch overhead + fills the GPU).
- [ ] Wire aggregation into the scan pipeline (both strands / all records per model in one launch).
- [ ] Forward kernel (last cascade stage; odds-space log-sum-exp, needs rescaling).
- [ ] Warp-per-window MSV — scoped to the ~10 M≥1000 models; or leave those on the CPU filter.
- [ ] Wire GPU Viterbi into the pipeline (after MSV survivors) once shared-mem Viterbi lands.
- [ ] Wire into `ornament-scfg`'s scan pipeline as an optional backend.
- [ ] Intra-DP (warp-per-window) for long models; multi-model batching.

## Feature gating

The `cuda` feature is **off by default**, so the crate — and the whole workspace — builds on
hosts with no CUDA toolchain (login nodes). Without it, the CPU oracle and batch builders are
available and `msv_nats_gpu` returns `GpuError::NotCompiled`. With `--features cuda`, `build.rs`
compiles `kernels/msv.cu` with `nvcc` and links the CUDA runtime.

## Building / testing the GPU path (this cluster)

CUDA 13.0 lives on the GPU nodes (`compgpu01-03`, 4× NVIDIA A30 = `sm_80`), not the login node.
Allocate a GPU with the `gpu_rbi` account on the `gpu` partition:

```bash
srun -A gpu_rbi -p gpu --gres=gpu:1 -c 8 -t 10 bash -lc '
  export CUDA_HOME=/usr/local/cuda-13.0
  export NVCC=$CUDA_HOME/bin/nvcc
  export PATH=$CUDA_HOME/bin:$PATH
  export CUDA_ARCH=sm_80        # A30 = Ampere
  export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
  cargo test -p ornament-gpu --features cuda -- --nocapture
'
```

`build.rs` env overrides: `NVCC` (path to nvcc), `CUDA_ARCH` (default `sm_80`), `CUDA_HOME`
(adds `$CUDA_HOME/{lib64,targets/x86_64-linux/lib}` to the link search path).

The parity test self-skips when no CUDA device is visible, so `--features cuda` is safe to run
anywhere; it only exercises the device when one is present.
