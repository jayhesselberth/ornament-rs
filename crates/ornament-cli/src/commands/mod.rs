//! CLI subcommand implementations.
//!
//! The engine layer (`search`, `scan`, `align`, `stat`) mirrors the Infernal
//! toolset (`cmsearch`/`cmscan`/`cmalign`/`cmstat`) over the native
//! `ornament-scfg` engine. The application layer (`trna`) is the
//! tRNA-modification analysis built on top.

pub mod align;
pub mod scan;
pub mod search;
pub mod stat;
pub mod trna;

use anyhow::{anyhow, Result};
use clap::ValueEnum;
use ornament_core::infernal::{write_tsv, CMHit};

/// Output format for hit reports (`search` / `scan`).
#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Format {
    /// Tab-separated hit table (the default; machine-friendly, pipeable).
    #[default]
    Tsv,
    /// Pretty-printed JSON array of hits.
    Json,
    /// Stockholm alignment of the hits (native engine only).
    Stockholm,
}

/// Search engine backing `search` / `scan`.
#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Engine {
    /// Native pure-Rust CM engine (no external binary or C toolchain). Default.
    #[default]
    Native,
    /// External `cmsearch` subprocess (Infernal must be on PATH). Oracle/fallback.
    Cmsearch,
}

/// E-value reporting threshold shared by the scan/search paths (mirrors `cmsearch -E`).
pub const E_VALUE: f64 = 1e-5;

/// Render a hit list as TSV or JSON (Stockholm is produced separately, from an MSA).
pub fn render_hits(format: Format, hits: &[CMHit]) -> Result<String> {
    match format {
        Format::Json => Ok(serde_json::to_string_pretty(hits)?),
        Format::Tsv => Ok(write_tsv(hits)),
        Format::Stockholm => Err(anyhow!(
            "stockholm output is assembled from alignments, not a hit list"
        )),
    }
}

/// Validate that a path exists, returning a clear error otherwise.
pub fn require_file(path: &str, what: &str) -> Result<()> {
    if !std::path::Path::new(path).exists() {
        return Err(anyhow!("{what} not found: {path}"));
    }
    Ok(())
}
