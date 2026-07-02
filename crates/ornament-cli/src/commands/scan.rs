//! `ornament scan` — scan many covariance models against sequences (the `cmscan` analog).
//!
//! This is the workhorse path: it runs every model in the `.cm` file (one model, or a whole
//! Rfam collection) over every FASTA record, parallelized across the model × record product.

use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use anyhow::{anyhow, Result};
use clap::Args;

use ornament_core::infernal::{
    scan_native_aligned, scan_native_multi, scan_native_multi_streaming, write_stockholm,
    write_tsv_rows, InfernalRunner, TSV_HEADER,
};

use super::{render_hits, require_file, resolve_model_path, Engine, Format, E_VALUE};
use crate::io::{emit, output_is_gz};
use crate::{progress, style};

#[derive(Args, Debug)]
pub struct ScanArgs {
    /// Input FASTA file
    #[arg(short, long)]
    pub input: String,

    /// Covariance model file (one model, or a whole Rfam collection)
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

pub fn run(args: ScanArgs) -> Result<()> {
    require_file(&args.input, "Input file")?;
    require_file(&args.cm, "CM file")?;

    // Prefer a pressed `<cm>.orm` sidecar for the native engine (skips per-run parse + banding); the
    // cmsearch subprocess always reads the `.cm` text directly.
    let cm = if args.engine == Engine::Native {
        resolve_model_path(&args.cm, !args.no_pressed)
    } else {
        args.cm.clone()
    };

    // Streaming fast path: native engine + TSV to stdout or a plain file. Hits are written as each
    // model finishes, so a long genome-vs-Rfam run shows progress, survives a crash, and never
    // holds the full hit set in memory. (`.gz` output needs a final gzip footer, so it falls
    // through to the buffered path below.)
    let can_stream = args.engine == Engine::Native
        && args.format == Format::Tsv
        && !output_is_gz(args.output.as_deref());
    if can_stream {
        return run_streaming(&args, &cm);
    }

    if args.format == Format::Stockholm {
        if args.engine != Engine::Native {
            return Err(anyhow!(
                "--format stockholm requires --engine native \
                 (the cmsearch subprocess does not produce native alignments)"
            ));
        }
        tracing::info!(
            "{} {} with {} (native engine, Stockholm alignment)",
            style::action("Scanning"),
            style::path(&args.input),
            style::path(&args.cm),
        );
        let spinner = progress::create_spinner("Scanning")?;
        let msas = scan_native_aligned(&cm, &args.input, E_VALUE)?;
        spinner.finish_and_clear();
        let n_hits: usize = msas.iter().map(|m| m.rows.len()).sum();
        tracing::info!(
            "Found {} hits across {} model(s)",
            style::count(n_hits),
            style::count(msas.len())
        );
        return emit(args.output.as_deref(), &write_stockholm(&msas));
    }

    let spinner = progress::create_spinner("Scanning")?;
    let hits = match args.engine {
        Engine::Native => {
            tracing::info!(
                "{} {} with {} (native engine, all models)",
                style::action("Scanning"),
                style::path(&args.input),
                style::path(&args.cm),
            );
            scan_native_multi(&cm, &args.input, E_VALUE)?
        }
        Engine::Cmsearch => {
            tracing::info!(
                "{} {} with {} (cmsearch subprocess)",
                style::action("Scanning"),
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

/// Native + TSV streaming path: write each model's hits as it completes. `cm` is the resolved model
/// path (a pressed `.orm` sidecar when one is current, else the `.cm`).
fn run_streaming(args: &ScanArgs, cm: &str) -> Result<()> {
    tracing::info!(
        "{} {} with {} (native engine, all models, streaming)",
        style::action("Scanning"),
        style::path(&args.input),
        style::path(&args.cm),
    );

    let mut w: Box<dyn Write + Send> = match args.output.as_deref() {
        Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    writeln!(w, "{TSV_HEADER}")?;
    let w = Mutex::new(w);

    let spinner = progress::create_spinner("Scanning")?;
    let total = AtomicUsize::new(0);

    let n = scan_native_multi_streaming(cm, &args.input, E_VALUE, |hits| {
        let mut s = String::new();
        write_tsv_rows(&mut s, hits);
        {
            let mut g = w.lock().unwrap();
            let _ = g.write_all(s.as_bytes());
            let _ = g.flush(); // flush per batch so partial results survive a crash
        }
        let seen = total.fetch_add(hits.len(), Ordering::Relaxed) + hits.len();
        spinner.set_message(format!("{seen} hits"));
    })?;

    spinner.finish_and_clear();
    w.into_inner().unwrap().flush()?;
    tracing::info!("Found {} hits", style::count(n));
    Ok(())
}
