# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ornament is a modification-aware tRNA scanner in Rust. It identifies tRNAs using Infernal covariance models, then analyzes sequence compatibility with known RNA modifications to flag "odd" tRNAs with modification-incompatible variants.

## Build Commands

The default build is pure Rust — it does **not** require `ext/` or a C toolchain.
The native port (`easel-rs`/`hmmer-rs`/`infernal-rs`) replaces the FFI; the legacy
`infernal-sys` FFI is optional behind `--features ffi` on `ornament-core`.

```bash
# Build (no ext/ needed)
cargo build

# Tests via nextest (preferred); falls back to `cargo test` if nextest absent
cargo nextest run --workspace
cargo test --doc            # nextest does not run doctests

# Run a specific crate's tests
cargo nextest run -p easel-rs

# Run the CLI
./target/debug/ornament mods --position 34
```

### Fast builds with mold + the differential-testing oracle (pixi)

The pixi `dev` environment provides the `mold` linker, `cargo-nextest`, and the
Infernal binaries (`cmsearch`/`cmstat`/`cmemit`) used as the differential-testing
oracle. Build/test wrappers prepend `mold -run` (LD_PRELOAD linker intercept):

```bash
pixi install -e dev              # one-time: solve mold + nextest + infernal
pixi run -e dev build            # mold -run cargo build
pixi run -e dev test             # mold -run cargo nextest run
pixi run -e dev check            # mold -run cargo check
pixi run -e dev clippy           # mold -run cargo clippy --workspace --all-targets
```

**Run builds/tests on SLURM, not the login node** — prefix with `srun -p rna -c 24`:

```bash
srun -p rna -c 24 pixi run -e dev test
srun -p rna -c 24 pixi run -e dev build
```

Only build `--release` when actually needed (release LTO is slow). The workspace
`default-members` excludes the legacy `infernal-sys` C-FFI crate, so bare `cargo`/`nextest`
commands skip it (build it explicitly with `-p infernal-sys` or `ornament-core`'s `ffi` feature).

`.cargo/config.toml` pins `target-cpu=x86-64-v3` (AVX2/FMA baseline, portable
across the cluster). AVX-512 hot kernels use runtime `is_x86_feature_detected!`
dispatch rather than a global baseline bump. See `.config/nextest.toml` for the
`default`/`ci` test profiles.

### Legacy C source (porting reference only)

```bash
# Clone Infernal/HMMER/Easel 1.1.5 into ext/ as a *reference* for the port
# (no longer linked into the default build).
./scripts/setup-deps.sh
```

## Architecture

### Workspace Structure

The project is migrating off the Infernal C FFI to a **native Rust port** of the
covariance-model search internals. Crates:

- **easel-rs** *(native)*: Rust port of the Easel subset — digital `Alphabet`
  (extended with first-class modified-base symbols), FASTA `Sequence` I/O,
  reverse-complement, and Gumbel/exponential survival stats for E-values.

- **hmmer-rs** *(native, scaffolding)*: Rust port of the HMMER p7 profile +
  MSV/Viterbi/Forward filters used as the CM acceleration filter (Phase 5).

- **infernal-rs** *(native, in progress)*: Rust port of the CM model, `.cm` parser,
  `cm_Configure`, emit-map, CYK/Inside scanning DP, alignment/traceback, tophits,
  and E-values. Emission arrays are sized from the alphabet's `K` so the model is
  alphabet-size-generic (modified-base discovery in Phase 8).

- **infernal-sys** *(legacy, optional)*: FFI bindings to the Infernal/HMMER/Easel C
  libraries. Compiled only with `--features ffi` on `ornament-core`; the `build.rs`
  compiles C from `ext/infernal/` via bindgen. The native crates replace this.

- **ornament-core**: Core library with modules for:
  - `modification/` - Types (`RnaBase`, `Modification`, `SprinzlPosition`), database of 12+ modifications with position expectations
  - `analysis/` - Compatibility scoring and odd tRNA detection
  - `infernal/` - Wrappers for CM search operations
  - `integration/` - modkit BedMethyl parsing
  - `output/` - JSON/TSV formatters

- **ornament-cli**: CLI binary with subcommands: `scan`, `analyze`, `compare`, `mods`

### Key Concepts

- **Sprinzl positions**: Standard tRNA numbering (1-76 with insertions like 17a, 20a). The `SprinzlMapper` maps CM alignment columns to these positions.

- **Modification compatibility**: Each modification has a `parent_base` (what it's derived from) and `incompatible_bases` (bases that preclude the modification). If a tRNA has an incompatible base at a position expecting a modification, it's flagged as "odd".

- **Modification-aware alphabet**: `easel-rs`'s `Alphabet` extends the canonical
  4-symbol RNA alphabet with registered modified bases (Ψ, D, I, m¹A, …). Each
  modified symbol's degeneracy vector points at its parent canonical base, so
  homology scoring is numerically identical to the stock alphabet while the
  modified identity is retained for annotation. This is the deliberate divergence
  from upstream Infernal.

- **Differential testing**: native DP/alignment output is validated against the
  `cmsearch` oracle (`cmsearch --max --nohmm` for exact unfiltered DP parity,
  `-A` for alignment parity). Oracle tests skip-gate when `cmsearch` is absent.

### External Dependencies (porting reference only)

The `ext/` directory (gitignored) holds the Infernal C source used as a **reference**
for the port — it is no longer linked into the default build:
- `ext/infernal/` - Infernal source (with `hmmer/` and `easel/` subdirectories)

Run `./scripts/setup-deps.sh --clean` to re-clone. Only needed for porting work or
the optional `--features ffi` legacy path.
