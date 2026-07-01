# ornament

A native, pure-Rust covariance-model (profile-SCFG) homology search engine — with a
tRNA-modification analysis application built on top.

Ornament is two layers:

1. **A native Rust reimplementation of Infernal.** A general-purpose covariance-model
   (profile-SCFG) homology search engine — *not* tRNA-specific. It reads standard `.cm`
   models and runs CYK/Inside scanning, alignment, and E-value calculation natively, with
   no C toolchain or FFI and no external `cmsearch`/`cmalign` binaries. It works with any
   Infernal covariance model (any Rfam family) and is validated for parity against the
   upstream `cmsearch` oracle. It is additionally **extended beyond Infernal** with a
   modification-aware alphabet (modified bases as first-class symbols).

2. **Ornament** — a tRNA-modification application. It identifies tRNAs with a covariance
   model, assigns Sprinzl positions via native alignment, and analyzes sequence
   compatibility with known RNA modifications to flag "odd" tRNAs that carry
   modification-incompatible variants.

## Features

- **Infernal-class CM toolkit**, native and pure-Rust — no `cmsearch` binary or C
  toolchain required:
  - `search` — one model vs. sequences (`cmsearch` analog)
  - `scan` — many models vs. sequences, streaming (`cmscan` analog)
  - `align` — align sequences to a model → Stockholm MSA (`cmalign` analog)
  - `stat` — model statistics (`cmstat` analog)
- **QDB-banded CYK/Inside** scanning DP, alignment/traceback, tophits, and Gumbel/
  exponential **E-values**
- **Profile-HMM (p7) filters** (MSV/Viterbi/Forward) for acceleration
- **Modification-aware alphabet** with modified bases (Ψ, D, I, m¹A, …) as first-class
  symbols
- **tRNA-modification analysis** (`trna` subcommand): Sprinzl mapping (1-76), a
  modification database (12+ modifications from MODOMICS), compatibility scoring, and
  modkit BedMethyl integration
- **Pipeable output**: TSV/JSON/Stockholm to stdout; status to stderr; gzip on `.gz`

## Installation

### Prerequisites

- Rust 1.96+ (with cargo)

That's it — the build is pure Rust. No C compiler, autoconf, `ext/`, or external Infernal
binaries are needed. (The `cmsearch` differential-testing oracle and the optional
`--engine cmsearch` subprocess path do require the Infernal toolchain on `PATH`; see
below.)

### Setup

```bash
git clone https://github.com/jayhesselberth/ornament-rs.git
cd ornament-rs
cargo build --release
```

The binary will be at `target/release/ornament`.

> The `ext/` Infernal/HMMER/Easel C sources (via `./scripts/setup-deps.sh`) are a
> **porting reference only** — they are not linked into the build.

## Usage

The engine commands mirror the Infernal toolset. TSV/JSON/Stockholm go to **stdout**
(plain, pipeable); status goes to **stderr**.

### Search one model against sequences (`cmsearch`)

```bash
# Native pure-Rust engine (default)
ornament search -c model.cm -i seqs.fa

# JSON to a file
ornament search -c model.cm -i seqs.fa -f json -o hits.json

# Alignment of the hits
ornament search -c model.cm -i seqs.fa -f stockholm
```

The `-f/--format` flag selects `tsv` (default), `json`, or `stockholm`. The
`--engine` flag selects the search backend:

| `--engine`  | Backend                                             | Requirements       |
|-------------|-----------------------------------------------------|--------------------|
| `native`    | pure-Rust `ornament-scfg` CYK scanner (**default**) | none               |
| `cmsearch`  | external `cmsearch` subprocess (oracle/fallback)    | Infernal on `PATH` |

### Scan many models against sequences (`cmscan`)

```bash
# Streams TSV as each model finishes (a whole Rfam collection is fine)
ornament scan -c Rfam.cm.gz -i genome.fa

# Gzip output (inferred from the .gz extension)
ornament scan -c Rfam.cm.gz -i genome.fa -o hits.tsv.gz

ornament scan -c Rfam.cm.gz -i genome.fa -f stockholm
```

### Align sequences to a model (`cmalign`)

```bash
ornament align -c model.cm -i members.fa            # Stockholm MSA to stdout
ornament align -c model.cm -i members.fa -o aln.sto
```

### Model statistics (`cmstat`)

```bash
ornament stat -c Rfam.cm.gz          # rounded table
ornament stat -c model.cm -f tsv     # pipeable TSV
```

### tRNA-modification analysis (`trna`)

The application layer lives under the `trna` subcommand.

```bash
# Show the modification database
ornament trna mods
ornament trna mods --position 34 --long

# Analyze modification compatibility of tRNA hits (JSON from search/scan)
ornament trna analyze --input trnas.json --threshold 0.8

# Compare with modkit modification calls
ornament trna compare --trna results.json --modkit mods.bedmethyl
```

> `analyze`/`compare`/`mods` also exist as hidden, deprecated top-level aliases that warn
> and dispatch into the `trna` namespace.

## How the tRNA analysis works

1. **Scan**: find tRNAs with the native `ornament-scfg` covariance-model engine (or, with
   `--engine cmsearch`, the external `cmsearch` subprocess).
2. **Map**: align hits and assign Sprinzl tRNA positions (1-76).
3. **Analyze**: check each position for modification compatibility — if a position expects
   modification X (derived from base Y) but the sequence carries a base Z incompatible with
   X, flag the tRNA as "odd".
4. **Report**: output compatibility scores and incompatibilities.

## Workspace structure

The engine is a native Rust port (no C FFI).

```
ornament-rs/
├── Cargo.toml                  # Workspace manifest
├── scripts/setup-deps.sh       # Clones ext/ Infernal C source (porting reference only)
├── ext/                        # Infernal C source — porting reference only (gitignored)
├── crates/
│   ├── ornament-alphabet/      # Modification-aware digital alphabet + FASTA I/O (noodles)
│   ├── ornament-stats/         # Gumbel/exponential E-value statistics
│   ├── ornament-hmm/           # Profile-HMM (p7) MSV/Viterbi/Forward filters
│   ├── ornament-scfg/          # Covariance-model (SCFG) engine — .cm parser, CYK/Inside
│   │                           #   scanning DP, alignment/traceback, tophits, E-values
│   ├── ornament-core/          # Core library
│   │   └── src/
│   │       ├── modification/   # Modification types and database
│   │       ├── analysis/       # Compatibility analysis
│   │       ├── infernal/       # Native + subprocess search wrappers
│   │       ├── integration/    # modkit BedMethyl parsing
│   │       └── output/         # JSON/TSV/Stockholm formatters
│   └── ornament-cli/           # `ornament` binary (commands/ + escapepod-style UX)
└── data/                       # Data files
```

## Modification database

The built-in database includes common eukaryotic tRNA modifications:

| Short Name | Full Name | Parent Base | Key Positions |
|------------|-----------|-------------|---------------|
| Psi | pseudouridine | U | 55 (universal) |
| D | dihydrouridine | U | 16, 17, 20 |
| m5U | 5-methyluridine | U | 54 |
| m1A | 1-methyladenosine | A | 58 |
| t6A | N6-threonylcarbamoyladenosine | A | 37 |
| i6A | N6-isopentenyladenosine | A | 37 |
| m1G | 1-methylguanosine | G | 37 |
| m7G | 7-methylguanosine | G | 46 |
| I | inosine | A | 34 (wobble) |
| Q | queuosine | G | 34 (wobble) |
| Cm | 2'-O-methylcytidine | C | 32 |
| m5C | 5-methylcytidine | C | 48 |

## Development

```bash
cargo build                 # debug build
cargo nextest run --workspace   # tests (falls back to `cargo test`)
cargo test --doc            # doctests (nextest does not run these)
cargo clippy
cargo fmt
```

Native DP/alignment output is validated against the `cmsearch` oracle
(`cmsearch --max --nohmm` for exact unfiltered DP parity, `-A` for alignment parity);
oracle tests skip-gate when `cmsearch` is absent. See `CLAUDE.md` for the pixi/mold/SLURM
build environment and criterion benchmarks.

## License

MIT

## References

- [Infernal](http://eddylab.org/infernal/) — RNA covariance model search
- [MODOMICS](https://genesilico.pl/modomics/) — RNA modification database
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) — tRNA detection
- [modkit](https://github.com/nanoporetech/modkit) — Nanopore modification toolkit
