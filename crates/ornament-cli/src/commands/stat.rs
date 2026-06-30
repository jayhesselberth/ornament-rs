//! `ornament stat` — report covariance-model statistics (the `cmstat` analog).

use anyhow::Result;
use clap::{Args, ValueEnum};
use tabled::builder::Builder;
use tabled::settings::Style;

use ornament_core::infernal::{cmstat, CmStat};

use super::require_file;
use crate::style;

/// Output format for `stat`.
#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum StatFormat {
    /// Human-readable rounded table (default).
    #[default]
    Table,
    /// Tab-separated values (one model per line; pipeable).
    Tsv,
    /// Pretty-printed JSON.
    Json,
}

#[derive(Args, Debug)]
pub struct StatArgs {
    /// Covariance model file (one model, or a whole collection)
    #[arg(short, long)]
    pub cm: String,

    /// Output format
    #[arg(short, long, value_enum, default_value_t = StatFormat::Table)]
    pub format: StatFormat,
}

const COLS: [&str; 11] = [
    "idx",
    "name",
    "accession",
    "clen",
    "states",
    "nodes",
    "W",
    "nseq",
    "GA",
    "TC",
    "calibrated",
];

fn opt_f32(v: Option<f32>) -> String {
    v.map_or_else(|| "-".to_string(), |x| format!("{x:.2}"))
}

fn row_cells(s: &CmStat) -> Vec<String> {
    vec![
        s.idx.to_string(),
        s.name.clone(),
        s.accession.clone().unwrap_or_else(|| "-".to_string()),
        s.clen.to_string(),
        s.states.to_string(),
        s.nodes.to_string(),
        s.w.to_string(),
        s.nseq.to_string(),
        opt_f32(s.ga),
        opt_f32(s.tc),
        if s.calibrated { "yes" } else { "no" }.to_string(),
    ]
}

pub fn run(args: StatArgs) -> Result<()> {
    require_file(&args.cm, "CM file")?;
    let stats = cmstat(&args.cm)?;

    match args.format {
        StatFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&stats)?);
        }
        StatFormat::Tsv => {
            println!("{}", COLS.join("\t"));
            for s in &stats {
                println!("{}", row_cells(s).join("\t"));
            }
        }
        StatFormat::Table => {
            eprintln!(
                "{} ({} model{})",
                style::header("Covariance model statistics"),
                stats.len(),
                if stats.len() == 1 { "" } else { "s" },
            );
            let mut builder = Builder::default();
            builder.push_record(COLS);
            for s in &stats {
                builder.push_record(row_cells(s));
            }
            println!("{}", builder.build().with(Style::rounded()));
        }
    }
    Ok(())
}
