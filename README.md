# ornament

Modification-aware tRNA scanner written in Rust.

Ornament identifies tRNAs in genomic sequences using Infernal covariance models, then analyzes sequence compatibility with known RNA modifications to flag "odd" tRNAs that have mutations incompatible with expected modifications.

## Features

- **tRNA scanning** via Infernal covariance models (FFI bindings to libinfernal)
- **Modification database** with 12+ tRNA modifications from MODOMICS
- **Sprinzl position mapping** (standard tRNA numbering 1-76)
- **Compatibility analysis** to detect modification-incompatible variants
- **modkit integration** for comparing with Nanopore modification calls

## Installation

### Prerequisites

- Rust 1.70+ (with cargo)
- C compiler (gcc or clang)
- autoconf

### Setup

1. Clone the repository:
```bash
git clone https://github.com/jayhesselberth/ornament.git
cd ornament
```

2. Set up dependencies (clones Infernal, HMMER, and Easel):
```bash
./scripts/setup-deps.sh
```

3. Build:
```bash
cargo build --release
```

The binary will be at `target/release/ornament`.

## Usage

### Show modification database

```bash
# List all modifications
ornament mods

# Show modifications expected at a specific position
ornament mods --position 34 --verbose

# Example output:
# Modifications expected at position 34:
#   inosine (I) - IsotypeSpecific conservation
#   queuosine (Q) - IsotypeSpecific conservation
```

### Scan sequences for tRNAs

```bash
ornament scan --input genome.fa --cm tRNA.cm --output results.json
```

### Analyze modification compatibility

```bash
ornament analyze --input trnas.json --threshold 0.8
```

### Compare with modkit calls

```bash
ornament compare --trna results.json --modkit mods.bedmethyl
```

## Project Structure

```
ornament/
├── Cargo.toml                 # Workspace manifest
├── scripts/
│   └── setup-deps.sh          # Dependency setup script
├── ext/                       # External deps (gitignored)
│   └── infernal/              # Infernal with hmmer/ and easel/
├── crates/
│   ├── infernal-sys/          # FFI bindings to Infernal/HMMER/Easel
│   │   ├── build.rs           # Compiles C libraries, generates bindings
│   │   └── src/lib.rs
│   ├── ornament-core/         # Core library
│   │   └── src/
│   │       ├── modification/  # Modification types and database
│   │       ├── analysis/      # Compatibility analysis
│   │       ├── infernal/      # Infernal interface
│   │       ├── integration/   # modkit integration
│   │       └── output/        # Output formatters
│   └── ornament-cli/          # CLI binary
│       └── src/main.rs
└── data/                      # Data files
```

## Modification Database

The database includes common eukaryotic tRNA modifications:

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

## How It Works

1. **Scan**: Use Infernal cmsearch to find tRNAs in input sequences
2. **Map**: Align hits to Sprinzl tRNA positions (1-76)
3. **Analyze**: Check each position for modification compatibility
   - If a position expects modification X (derived from base Y)
   - But the sequence has base Z (incompatible with X)
   - Flag as "odd" tRNA
4. **Report**: Output compatibility scores and incompatibilities

## Development

```bash
# Run tests
cargo test

# Build in debug mode
cargo build

# Run clippy
cargo clippy

# Format code
cargo fmt
```

## License

MIT

## References

- [Infernal](http://eddylab.org/infernal/) - RNA covariance model search
- [MODOMICS](https://genesilico.pl/modomics/) - RNA modification database
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) - tRNA detection
- [modkit](https://github.com/nanoporetech/modkit) - Nanopore modification toolkit
