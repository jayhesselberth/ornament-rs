//! ornament — a native covariance-model (profile-SCFG) homology search engine,
//! with a tRNA-modification analysis application.
//!
//! The engine layer (`search`, `scan`, `align`, `stat`) mirrors the Infernal
//! toolset over the pure-Rust `ornament-scfg` engine; the `trna` layer is the
//! modification-analysis application built on top.

use anyhow::Result;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use clap::{Parser, Subcommand};
use tracing_subscriber::EnvFilter;

mod commands;
mod io;
mod progress;
mod style;

use commands::{align, press, scan, search, stat, trna};

/// clap help palette (mirrors the runtime status palette): bold-yellow headers/usage,
/// bold-green literals, cyan placeholders.
const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Yellow.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default());

#[derive(Parser)]
#[command(name = "ornament")]
#[command(author, version, styles = STYLES)]
#[command(
    about = "A native covariance-model homology search engine, with tRNA-modification analysis"
)]
#[command(
    long_about = "ornament is a native, pure-Rust covariance-model (profile-SCFG) homology \
search engine — a reimplementation of the Infernal toolset that needs no external binary or C \
toolchain. It works with any Infernal covariance model (any Rfam family).\n\n\
The engine commands mirror Infernal: `search` (cmsearch), `scan` (cmscan), `align` (cmalign), and \
`stat` (cmstat). The `trna` subcommand is a tRNA-modification analysis application built on the \
engine."
)]
#[command(after_help = "\
Examples:
  ornament search -c model.cm -i seqs.fa          Search one model against sequences
  ornament scan   -c Rfam.cm.gz -i genome.fa      Scan a model collection (cmscan-style)
  ornament align  -c model.cm -i members.fa       Align sequences to a model (Stockholm)
  ornament stat   -c Rfam.cm.gz                    Show model statistics
  ornament trna analyze -i hits.json              Analyze tRNA modification compatibility
  ornament trna mods --position 34                Show modifications expected at a position

Engine commands write TSV/JSON/Stockholm to stdout (pipeable); status goes to stderr.
Color follows NO_COLOR / CLICOLOR / TTY detection.
")]
struct Cli {
    /// Silence all output except errors. Overrides `--verbose`.
    #[arg(short = 'q', long, global = true)]
    quiet: bool,

    /// Increase log verbosity. Status messages show by default (info);
    /// `-v` = debug, `-vv` = trace. `RUST_LOG` takes precedence when set.
    #[arg(short = 'v', long, global = true, action = clap::ArgAction::Count)]
    verbose: u8,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Search one covariance model against sequences (cmsearch analog)
    #[command(after_help = "\
Examples:
  ornament search -c model.cm -i seqs.fa
  ornament search -c model.cm -i seqs.fa -f json -o hits.json
  ornament search -c model.cm -i seqs.fa -f stockholm     Alignment of the hits")]
    Search(search::SearchArgs),

    /// Scan many covariance models against sequences (cmscan analog)
    #[command(after_help = "\
Examples:
  ornament scan -c Rfam.cm.gz -i genome.fa                Streams TSV as each model finishes
  ornament scan -c Rfam.cm.gz -i genome.fa -o hits.tsv.gz Gzip output
  ornament scan -c Rfam.cm.gz -i genome.fa -f stockholm")]
    Scan(scan::ScanArgs),

    /// Align given sequences to a model, emit a Stockholm MSA (cmalign analog)
    #[command(after_help = "\
Examples:
  ornament align -c model.cm -i members.fa
  ornament align -c model.cm -i members.fa -o aln.sto")]
    Align(align::AlignArgs),

    /// Report covariance-model statistics (cmstat analog)
    #[command(after_help = "\
Examples:
  ornament stat -c Rfam.cm.gz
  ornament stat -c model.cm -f tsv")]
    Stat(stat::StatArgs),

    /// Prebuild a covariance-model database, precomputing bands (cmpress analog)
    #[command(after_help = "\
Examples:
  ornament press -c Rfam.cm.gz            Writes Rfam.cm.gz.orm alongside the model
  ornament press -c model.cm -f           Overwrite an existing pressed database

`scan`/`search` auto-detect the `<cm>.orm` sidecar and skip re-parsing + re-banding each run.
Pass `--no-pressed` to those commands to ignore a sidecar.")]
    Press(press::PressArgs),

    /// tRNA-modification analysis application (analyze, compare, mods)
    Trna {
        #[command(subcommand)]
        command: trna::TrnaCommand,
    },

    /// Deprecated: use `ornament trna analyze`
    #[command(hide = true)]
    Analyze(trna::AnalyzeArgs),
    /// Deprecated: use `ornament trna compare`
    #[command(hide = true)]
    Compare(trna::CompareArgs),
    /// Deprecated: use `ornament trna mods`
    #[command(hide = true)]
    Mods(trna::ModsArgs),
}

/// Install the tracing subscriber: terse status lines on stderr, level/verbosity gated by
/// `-q`/`-v` (and `RUST_LOG`), color following the shared `NO_COLOR`/`CLICOLOR`/TTY gate.
fn init_logging(quiet: bool, verbose: u8) {
    let default = match (quiet, verbose) {
        (true, _) => "error",
        (_, 0) => "info",
        (_, 1) => "debug",
        _ => "trace",
    };
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(default));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_ansi(style::use_color())
        .without_time()
        .with_target(false)
        .init();
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.quiet, cli.verbose);

    match cli.command {
        Commands::Search(args) => search::run(args),
        Commands::Scan(args) => scan::run(args),
        Commands::Align(args) => align::run(args),
        Commands::Stat(args) => stat::run(args),
        Commands::Press(args) => press::run(args),
        Commands::Trna { command } => trna::run(command),

        // Deprecated top-level aliases — dispatch to the `trna` namespace with a notice.
        Commands::Analyze(args) => {
            tracing::warn!("`ornament analyze` is deprecated; use `ornament trna analyze`");
            trna::run(trna::TrnaCommand::Analyze(args))
        }
        Commands::Compare(args) => {
            tracing::warn!("`ornament compare` is deprecated; use `ornament trna compare`");
            trna::run(trna::TrnaCommand::Compare(args))
        }
        Commands::Mods(args) => {
            tracing::warn!("`ornament mods` is deprecated; use `ornament trna mods`");
            trna::run(trna::TrnaCommand::Mods(args))
        }
    }
}
