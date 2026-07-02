//! `ornament search` — search a single covariance model against sequences (the `cmsearch` analog).

use anyhow::{anyhow, Result};
use clap::Args;

use ornament_core::infernal::{scan_native, scan_native_aligned, write_stockholm, InfernalRunner};

use super::{render_hits, require_file, resolve_model_path, Engine, Format, E_VALUE};
use crate::io::emit;
use crate::{progress, style};

#[derive(Args, Debug)]
pub struct SearchArgs {
    /// Input FASTA file
    #[arg(short, long)]
    pub input: String,

    /// Covariance model file (single model)
    #[arg(short, long)]
    pub cm: String,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<String>,

    /// Output format
    #[arg(short, long, value_enum, default_value_t = Format::Tsv)]
    pub format: Format,

    /// Search engine: native (pure Rust, default) or cmsearch (subprocess)
    #[arg(long, value_enum, default_value_t = Engine::Native)]
    pub engine: Engine,

    /// Ignore any pressed `<cm>.orm` sidecar and parse the `.cm` directly
    #[arg(long)]
    pub no_pressed: bool,
}

pub fn run(args: SearchArgs) -> Result<()> {
    require_file(&args.input, "Input file")?;
    require_file(&args.cm, "CM file")?;

    // Prefer a pressed `<cm>.orm` sidecar for the native engine; cmsearch reads the `.cm` directly.
    let cm = if args.engine == Engine::Native {
        resolve_model_path(&args.cm, !args.no_pressed)
    } else {
        args.cm.clone()
    };

    if args.format == Format::Stockholm {
        if args.engine != Engine::Native {
            return Err(anyhow!(
                "--format stockholm requires --engine native \
                 (the cmsearch subprocess does not produce native alignments)"
            ));
        }
        tracing::info!(
            "{} {} with {} (native engine, Stockholm alignment)",
            style::action("Searching"),
            style::path(&args.input),
            style::path(&args.cm),
        );
        let msas = scan_native_aligned(&cm, &args.input, E_VALUE)?;
        let n_hits: usize = msas.iter().map(|m| m.rows.len()).sum();
        tracing::info!("Found {} hits", style::count(n_hits));
        return emit(args.output.as_deref(), &write_stockholm(&msas));
    }

    let spinner = progress::create_spinner("Searching")?;
    let hits = match args.engine {
        Engine::Native => {
            tracing::info!(
                "{} {} with {} (native engine)",
                style::action("Searching"),
                style::path(&args.input),
                style::path(&args.cm),
            );
            scan_native(&cm, &args.input, E_VALUE)?
        }
        Engine::Cmsearch => {
            tracing::info!(
                "{} {} with {} (cmsearch subprocess)",
                style::action("Searching"),
                style::path(&args.input),
                style::path(&args.cm),
            );
            InfernalRunner::new()
                .with_cm(&args.cm)
                .with_e_value(E_VALUE)
                .cmsearch(&args.input)?
        }
    };
    spinner.finish_and_clear();

    tracing::info!("Found {} hits", style::count(hits.len()));
    emit(args.output.as_deref(), &render_hits(args.format, &hits)?)
}
