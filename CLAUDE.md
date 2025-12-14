# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ornament is a modification-aware tRNA scanner in Rust. It identifies tRNAs using Infernal covariance models, then analyzes sequence compatibility with known RNA modifications to flag "odd" tRNAs with modification-incompatible variants.

## Build Commands

```bash
# First-time setup (clones Infernal/HMMER/Easel to ext/)
./scripts/setup-deps.sh

# Build
cargo build

# Run tests
cargo test

# Run a specific test
cargo test test_name

# Run tests in a specific crate
cargo test -p ornament-core

# Run the CLI
./target/debug/ornament mods --position 34
```

## Architecture

### Workspace Structure

Three crates in a Cargo workspace:

- **infernal-sys**: FFI bindings to Infernal/HMMER/Easel C libraries. The `build.rs` compiles the C code from `ext/infernal/` and generates Rust bindings via bindgen. Statically links `libinfernal.a`, `libhmmer.a`, `libeasel.a`.

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

- **FFI build**: The infernal-sys build.rs runs `./configure && make` in the build directory. Uses `infernal-1.1.5` tags for synchronized Infernal/HMMER/Easel versions.

### External Dependencies

The `ext/` directory (gitignored) contains:
- `ext/infernal/` - Infernal source (with `hmmer/` and `easel/` subdirectories)

Run `./scripts/setup-deps.sh --clean` to re-clone dependencies.
