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
- [ ] Warp-per-window for long models — the 103 fallback models are the measured bottleneck.
- [ ] Reusable device context — hoist stream/pinned-buffer setup out of the per-call path so
      many-model scans don't pay allocation per model.
- [ ] Viterbi + Forward kernels (rest of the cascade).
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
