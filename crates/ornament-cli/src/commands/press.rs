//! `ornament press` — prebuild a covariance-model database (the `cmpress` analog).
//!
//! Prepares every model in a `.cm(.gz)` once — configure, QDB bands, filter HMM — and writes a
//! binary `<cm>.orm` sidecar. Subsequent `scan`/`search` runs auto-detect the sidecar and skip that
//! per-run preparation (parse + configure + band computation), which is the bulk of a many-model
//! scan's startup cost. Most useful when the same model collection (e.g. Rfam) is searched repeatedly.

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;

use ornament_core::infernal::press_cm;

use super::require_file;
use crate::style;

#[derive(Args, Debug)]
pub struct PressArgs {
    /// Covariance model file to press (one model, or a whole Rfam collection)
    #[arg(short, long)]
    pub cm: String,

    /// Output path for the pressed database (default: `<cm>.orm`; `.gz` suffix compresses it)
    #[arg(short, long)]
    pub output: Option<String>,

    /// Overwrite an existing pressed database
    #[arg(short, long)]
    pub force: bool,
}

pub fn run(args: PressArgs) -> Result<()> {
    require_file(&args.cm, "CM file")?;

    let out: Option<PathBuf> = args.output.as_deref().map(PathBuf::from);
    tracing::info!("{} {}", style::action("Pressing"), style::path(&args.cm),);

    let spinner = crate::progress::create_spinner("Pressing")?;
    let (n, written) = press_cm(&args.cm, out.as_deref(), args.force)?;
    spinner.finish_and_clear();

    tracing::info!(
        "Pressed {} model(s) → {}",
        style::count(n),
        style::path(written.display()),
    );
    Ok(())
}
