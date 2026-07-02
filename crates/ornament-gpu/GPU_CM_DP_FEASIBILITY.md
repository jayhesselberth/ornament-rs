# GPU CM covariance-model DP — feasibility scoping

Motivation: the realistic all-Rfam-vs-genome scan is ~99% CM DP, not p7 filters (setup ≈ 3 s of
354 s; GPU MSV only bought ~3%). So the only lever that meaningfully speeds the scan is the CM
covariance-model DP itself — the profile-SCFG CYK/Inside scan (`search.rs::scan_core`) and the
per-hit alignment (`align.rs`). This document scopes whether that DP can go on the GPU.

## What the DP actually is (from `search.rs::scan_core`)

A banded profile-SCFG DP over three axes: **state `v`**, **end position `j`**, **subsequence
length `d`** — `alpha[v][j][d]`. Semiring is max-plus (CYK) or log-sum-exp (Inside).

Loop structure (the shipping, optimized version):
```
for j in 1..=L:                      # SEQUENTIAL (j reads j-1 and same-j cells)
  for v in (M-1)..=1:                # reverse topological order (child before parent)
    for d in [dmin[v] .. dmax[v]]:   # QDB-banded, RAGGED per state
      alpha[v][j][d] = ⊕_children (alpha[child][.][d-sd] ⊗ tsc) ⊗ emission
```
Cell kinds:
- **Emitting / D / S states**: a child-transition fold over `d` (already SIMD-vectorized on CPU via
  `S::accumulate`) + a per-`d` emission. Regular, parallel over `d`.
- **Bifurcation (B) states**: `alpha[B][j][d] = ⊕_k left(j-k, d-k) ⊗ right(j, k)` — an O(d)
  convolution using a `W+1` rolling "BEGL deck". Irregular; variable count per model.
- **IL insert self-loops**: `alpha[v][j][d]` reads `alpha[v][j][d-1]` — **loop-carried in `d`**,
  cannot be vectorized; a sequential recurrence.
- A **gamma semi-HMM** pass + traceback resolves non-overlapping hits (sequential, pointer-chasing).

Cost ≈ O(L · Σ_v band_width_v) for regular states, plus O(L · Σ_B band² ) for bifurcations. For a
large structured RNA (SSU rRNA: M≈1500, W≈1500) this is the scan's dominant cost.

## Parallelism map (what fits the GPU, what doesn't)

| axis | parallel? | notes |
|------|-----------|-------|
| `d` (width) | **yes** | the natural thread/SIMD axis; already vectorized on CPU |
| `v` (states) | **partial** | reverse-topo dependency (child→parent); only *sibling subtrees* are independent — extracting that needs a topological-level decomposition of the guide tree |
| `j` (position) | **no** | fundamentally sequential: `j` reads `j-1` (prv row) and same-`j` left-emitter cells |
| across windows/models | **yes** | many independent CM scans could batch (like the filter batching) |

So a single CM scan has intra-`j` parallelism (over `d`, and partially over `v`) but a **sequential
outer `j` loop of length W** — the same shape as GPU RNA folding (Nussinov/Zuker), which is known
to be **hard to make efficient** because the sequential dimension throttles occupancy and each
step's parallel work shrinks with the banding.

## Memory

Rolling 2-row `alpha` = `2·M·(W+1)` floats (≈18 MB at M=W=1500) + the BEGL deck `O(W · Σ_B dmax)`.
The pipeline already skips models whose DP would exceed 1.5 GiB, so the *feasible* models fit in an
A30's 24 GB. Memory is **not** the blocker (banded); the giant rRNA models that don't fit are
already CPU-skipped.

## The hard parts (ranked)

1. **Sequential `j` (length W)** — caps GPU utilization; each step's parallel work is only
   `Σ_v band_v`, which shrinks under banding. This is intrinsic to CYK and the main efficiency limit.
2. **State-dependency wavefront** — within a `j`-step, states aren't all independent; naive impl
   serializes M states. Real speedup needs guide-tree level decomposition (independent subtrees per
   wavefront level) — significant complexity.
3. **Bifurcation convolutions** — O(d) reductions with the rolling deck; irregular, model-dependent.
4. **IL self-loops** — loop-carried in `d`, inherently sequential (a scan/prefix pattern).
5. **Ragged QDB bands** — per-state `[dmin,dmax]` ⇒ warp divergence + irregular memory; the exact
   feature that makes the CPU fast makes the GPU awkward.
6. **Gamma resolution + traceback + per-hit alignment** — sequential, pointer-chasing; stays on CPU
   (so every scan still round-trips CPU↔GPU).

## Precedent

GPU SCFG/CYK parsing and GPU RNA folding (CUDA Nussinov/Zuker, a few research GPU-Infernal efforts)
exist and show *modest* speedups with *large* effort. Notably, **Infernal ships CPU-only** for the
CM DP; the GPU attention in the field went to the HMMER p7 filters (which we've now done) precisely
because the linear filter DPs are easy and the SCFG DP is not. There is **no production GPU
cmsearch**.

## Verdict

**Technically feasible, but a poor bet:** high effort (weeks — a full banded CYK+Inside port with
bifurcations, IL self-loops, ragged bands, and parity validation against the CPU oracle), high risk
(the sequential-`j` + wavefront + banding structure limits achievable speedup), and a payoff bounded
to the few heavy models per scan. The effort/reward is far worse than the p7 filters were.

### Better-scoped options, in decreasing attractiveness

1. **CPU-side heavy-model optimization (recommended).** The scan cost is a *few* large models
   (rRNA). Tighter scan banding (a larger `beta` = looser tail loss ⇒ narrower bands, for the *scan*
   only, keeping the calibrated align beta), or capping/streaming the rRNA DP, would cut the actual
   bottleneck with far less risk. The codebase already has the memory-guard + banded-align machinery
   to build on. This is the highest reward-per-effort path to a faster scan.
2. **Batch many *small*-model CM scans on GPU** (one block per (model, survivor-window) DP), reusing
   the filter-batching + tile-aggregation infrastructure already built. But small models are already
   cheap on CPU, so this speeds the part that isn't the bottleneck — low value.
3. **Full GPU banded CYK for large models** — the direct but hardest path; only worth it as a
   research effort, and even then the sequential `j` likely caps the win below what the effort costs.

**Recommendation:** do **not** invest in GPU CM-DP now. The realistic-scan measurement says the
p7-filter GPU work (done, validated) is the sensible extent of GPU offload for this engine; further
scan speedup should come from CPU-side heavy-model DP optimization, which targets the true
bottleneck at a fraction of the risk.
