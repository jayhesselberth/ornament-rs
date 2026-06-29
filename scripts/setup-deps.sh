#!/bin/bash
# Clone the Infernal/HMMER/Easel 1.1.5 C sources into ext/.
#
# PORTING REFERENCE ONLY: ext/ is *not* a build dependency. The default build is
# pure Rust (easel-rs/hmmer-rs/infernal-rs) and needs neither ext/ nor a C
# toolchain. These sources exist only as a reference for the native port and for
# the optional legacy C-FFI path (ornament-core --features ffi / -p infernal-sys).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
EXT_DIR="$ROOT_DIR/ext"

# Default to stable infernal-1.1.5 release (synchronized across all repos)
# Use develop branches only when all three are compatible
INFERNAL_BRANCH="${INFERNAL_BRANCH:-infernal-1.1.5}"
HMMER_BRANCH="${HMMER_BRANCH:-infernal-1.1.5}"
EASEL_BRANCH="${EASEL_BRANCH:-infernal-1.1.5}"

# Check for --clean flag to remove existing dependencies
if [ "$1" = "--clean" ]; then
    echo "Cleaning existing dependencies..."
    rm -rf "$EXT_DIR/infernal"
fi

echo "Setting up ornament dependencies..."
echo "  Infernal branch: $INFERNAL_BRANCH"
echo "  HMMER branch: $HMMER_BRANCH"
echo "  Easel branch: $EASEL_BRANCH"

mkdir -p "$EXT_DIR"

# Clone Infernal
if [ ! -d "$EXT_DIR/infernal" ]; then
    echo "Cloning Infernal..."
    git clone -b "$INFERNAL_BRANCH" https://github.com/EddyRivasLab/infernal.git "$EXT_DIR/infernal"
else
    echo "Infernal already exists, skipping..."
fi

# Clone HMMER into infernal directory
if [ ! -d "$EXT_DIR/infernal/hmmer" ]; then
    echo "Cloning HMMER..."
    git clone -b "$HMMER_BRANCH" https://github.com/EddyRivasLab/hmmer.git "$EXT_DIR/infernal/hmmer"
else
    echo "HMMER already exists, skipping..."
fi

# Clone Easel into infernal directory
if [ ! -d "$EXT_DIR/infernal/easel" ]; then
    echo "Cloning Easel..."
    git clone -b "$EASEL_BRANCH" https://github.com/EddyRivasLab/easel.git "$EXT_DIR/infernal/easel"
else
    echo "Easel already exists, skipping..."
fi

# Run autoconf if configure doesn't exist
if [ ! -f "$EXT_DIR/infernal/configure" ]; then
    echo "Running autoconf..."
    cd "$EXT_DIR/infernal"
    autoconf
fi

echo ""
echo "Porting-reference sources cloned into ext/."
echo "Note: the default 'cargo build' is pure Rust and does NOT use ext/."
echo "ext/ is only needed for porting work or the optional 'ffi' feature."
