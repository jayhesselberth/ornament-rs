# ornament-gpu

CUDA-accelerated p7 filter kernels — a **prototype** toward offloading the all-model /
many-genomes `scan` to GPU.

The scan's heaviest, most parallel workload is the p7 prefilter cascade run over a huge batch
of `(model, window)` tiles. This crate begins with the cheapest, most GPU-friendly stage — the
**MSV filter** — batched **one thread per window** on the device (`kernels/msv.cu`). The kernel
is a faithful f32 max-plus transcription of the scalar [`ornament_hmm::msv_nats`] oracle, so the
device result matches the CPU to f32 rounding (verified by the `gpu_matches_cpu` parity test).

## Status

- [x] Batch MSV kernel (f32, one-thread-per-window) + CUDA-Runtime-API host wrappers.
- [x] Feature-gated build (`nvcc` via `build.rs`), CPU oracle, GPU↔CPU parity test.
- [ ] uint8 reduced-precision MSV (throughput; mirrors `ornament-hmm::msv_simd`).
- [ ] Viterbi + Forward kernels (rest of the cascade).
- [ ] Wire into `ornament-scfg`'s scan pipeline as an optional backend.
- [ ] Intra-DP parallelism / shared-memory staging for long models; multi-model batching.

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
