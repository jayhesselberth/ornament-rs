//! CLI subcommand implementations.
//!
//! The engine layer (`search`, `scan`, `align`, `stat`) mirrors the Infernal
//! toolset (`cmsearch`/`cmscan`/`cmalign`/`cmstat`) over the native
//! `ornament-scfg` engine. The application layer (`trna`) is the
//! tRNA-modification analysis built on top.

pub mod align;
pub mod press;
pub mod scan;
pub mod search;
pub mod stat;
pub mod trna;

use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};
use clap::ValueEnum;
use ornament_core::infernal::{write_tsv, CMHit};

/// Pressed-DB sidecar extension (`Rfam.cm` → `Rfam.cm.orm`). Mirrors `ornament_scfg::PRESS_EXT`.
const PRESS_EXT: &str = "orm";

/// Resolve which model file a native `scan`/`search` should actually open. When `use_pressed` is
/// set (the default) and a pressed sidecar `<cm>.orm` exists and is at least as new as the `.cm`
/// it was built from, the sidecar is returned so the run skips re-parsing + re-banding; a *stale*
/// sidecar (older than the `.cm`) is ignored with a warning. A pressed file passed directly is used
/// as-is. `use_pressed = false` (`--no-pressed`) always returns `cm` unchanged.
pub fn resolve_model_path(cm: &str, use_pressed: bool) -> String {
    if !use_pressed {
        return cm.to_string();
    }
    let cm_path = Path::new(cm);
    // The user passed a pressed file directly — the core loader detects it by content.
    if cm_path.extension().and_then(|e| e.to_str()) == Some(PRESS_EXT) {
        return cm.to_string();
    }
    let sidecar = press_sidecar(cm_path);
    if !sidecar.exists() {
        return cm.to_string();
    }
    if sidecar_is_stale(cm_path, &sidecar) {
        tracing::warn!(
            "pressed DB {} is older than {}; ignoring it (re-run `ornament press`)",
            sidecar.display(),
            cm
        );
        return cm.to_string();
    }
    tracing::info!("using pressed CM database {}", sidecar.display());
    sidecar.to_string_lossy().into_owned()
}

/// The default sidecar path for a `.cm`: append `.orm` to the full name.
fn press_sidecar(cm: &Path) -> PathBuf {
    let mut s = cm.as_os_str().to_os_string();
    s.push(".");
    s.push(PRESS_EXT);
    PathBuf::from(s)
}

/// Whether the pressed sidecar predates the `.cm` it was built from (so it may be out of date).
/// If either mtime can't be read, assume fresh (don't cry wolf).
fn sidecar_is_stale(cm: &Path, sidecar: &Path) -> bool {
    match (
        cm.metadata().and_then(|m| m.modified()),
        sidecar.metadata().and_then(|m| m.modified()),
    ) {
        (Ok(cm_t), Ok(side_t)) => side_t < cm_t,
        _ => false,
    }
}

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
