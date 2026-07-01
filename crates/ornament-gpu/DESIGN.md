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
- **Overlap:** CUDA streams + double-buffered pinned host buffers. While the GPU scores batch *i*,
  the host builds/copies batch *i+1*. Target: keep the A30 saturated, never blocked on H2D/D2H.

### Level 2 — on-device funnel (stream compaction)
- MSV kernel writes a survivor mask → compact → Viterbi kernel consumes only survivors → compact
  → Forward. Compaction between kernels keeps later, costlier kernels dense (no warps full of
  already-pruned lanes). Survivors flow back to the host for CYK.

## Concrete wins, in priority order

1. **Resident sequence + offset windows.** Upload each strand to the device *once*; make windows
   `(start,end)` offsets into device-resident memory instead of copying residues per window. The
   scan tiles at length `2W` step `W` (2× overlap), so the current per-window copy moves each base
   ~2×. This likely beats the uint8 conversion on its own.
2. **uint8 MSV** *(done — `msv_u8_batch_kernel`)*: 4× smaller emission table → better bandwidth on
   the memory-bound one-thread-per-window batch.
3. **Batch across models and sequences.** A single small model (tRNA, M≈71) underutilizes the GPU.
   - *Many-genomes, one model:* one resident emission table + a huge window batch. Current design.
   - *All-model, one genome (~4000 Rfam):* iterate models but keep emission tables resident/cached
     and pack small models together (per-window model-id index into resident tables).
4. **Kernel shape by model length.** One-thread-per-window is ideal for short models but spills
   scratch for long M. Endgame: one-**warp**-per-window with a shared-memory DP row for large M;
   dispatch on M.

## Parity discipline

Every device kernel is validated against a trusted CPU oracle to f32 rounding:
- `msv_batch_kernel` (f32) ↔ `ornament_hmm::msv_nats`.
- `msv_u8_batch_kernel` (uint8) ↔ `ornament_hmm::MsvProfile::score_bits` (the striped uint8 filter).

The CPU path stays the reference; the GPU only changes speed, never the reported result.
