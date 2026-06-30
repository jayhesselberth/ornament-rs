//! `ornament align` — align given sequences to a model and emit a Stockholm MSA (the `cmalign`
//! analog). Performs no search: every input sequence is aligned to the model over its full length.

use anyhow::Result;
use clap::Args;

use ornament_core::infernal::{align_sequences, write_stockholm};

use super::require_file;
use crate::io::emit;
use crate::{progress, style};

#[derive(Args, Debug)]
pub struct AlignArgs {
    /// Input FASTA file (sequences to align — assumed members of the model)
    #[arg(short, long)]
    pub input: String,

    /// Covariance model file (single model)
    #[arg(short, long)]
    pub cm: String,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<String>,
}

pub fn run(args: AlignArgs) -> Result<()> {
    require_file(&args.input, "Input file")?;
    require_file(&args.cm, "CM file")?;

    tracing::info!(
        "{} {} to {} (native engine)",
        style::action("Aligning"),
        style::path(&args.input),
        style::path(&args.cm),
    );

    let spinner = progress::create_spinner("Aligning")?;
    let msa = align_sequences(&args.cm, &args.input)?;
    spinner.finish_and_clear();

    tracing::info!("Aligned {} sequences", style::count(msa.rows.len()));
    emit(
        args.output.as_deref(),
        &write_stockholm(std::slice::from_ref(&msa)),
    )
}
