//! ornament CLI - Modification-aware tRNA scanner
//!
//! Scans genomic sequences for tRNAs and analyzes modification compatibility.

use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Parser)]
#[command(name = "ornament")]
#[command(author, version, about = "Modification-aware tRNA scanner", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Scan sequences for tRNAs
    Scan {
        /// Input FASTA file
        #[arg(short, long)]
        input: String,

        /// Covariance model file
        #[arg(short, long)]
        cm: Option<String>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<String>,

        /// Output format (json, tsv)
        #[arg(short, long, default_value = "json")]
        format: String,
    },

    /// Analyze modification compatibility of tRNA sequences
    Analyze {
        /// Input tRNA sequences (FASTA or JSON from scan)
        #[arg(short, long)]
        input: String,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<String>,

        /// Compatibility score threshold for "odd" tRNAs
        #[arg(short, long, default_value = "0.8")]
        threshold: f64,
    },

    /// Compare with modkit modification calls
    Compare {
        /// tRNA analysis results (JSON)
        #[arg(short, long)]
        trna: String,

        /// modkit BedMethyl file
        #[arg(short, long)]
        modkit: String,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<String>,
    },

    /// Show modification database
    Mods {
        /// Filter by position (e.g., "34", "55")
        #[arg(short, long)]
        position: Option<String>,

        /// Show detailed information
        #[arg(short, long)]
        verbose: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Scan { input, cm, output: _, format: _ } => {
            println!("Scanning {} for tRNAs...", input);
            if let Some(cm_path) = cm {
                println!("Using CM: {}", cm_path);
            }
            // TODO: Implement scanning
            println!("Scan not yet implemented");
        }

        Commands::Analyze { input, output: _, threshold } => {
            println!("Analyzing modification compatibility in {}...", input);
            println!("Threshold: {}", threshold);
            // TODO: Implement analysis
            println!("Analysis not yet implemented");
        }

        Commands::Compare { trna, modkit, output: _ } => {
            println!("Comparing {} with modkit calls from {}...", trna, modkit);
            // TODO: Implement comparison
            println!("Comparison not yet implemented");
        }

        Commands::Mods { position, verbose } => {
            let db = ornament_core::modification::ModificationDatabase::eukaryotic();

            if let Some(pos) = position {
                let sprinzl = ornament_core::SprinzlPosition(pos.clone());
                let expectations = db.get_expectations(&sprinzl);
                if !expectations.is_empty() {
                    println!("Modifications expected at position {}:", pos);
                    for exp in expectations {
                        for modification in &exp.modifications {
                            if verbose {
                                println!("  {} ({}) - {:?} conservation",
                                         modification.name,
                                         modification.short_name,
                                         exp.conservation);
                            } else {
                                println!("  {}", modification.short_name);
                            }
                        }
                    }
                } else {
                    println!("No modifications expected at position {}", pos);
                }
            } else {
                println!("Modification database ({} modifications):", db.modifications().len());
                for (name, modification) in db.modifications() {
                    if verbose {
                        println!("  {} ({}) - derived from {:?}, genomic: {:?}",
                                 name,
                                 modification.short_name,
                                 modification.parent_base,
                                 modification.genomic_expectation);
                    } else {
                        println!("  {} ({})", modification.short_name, name);
                    }
                }
            }
        }
    }

    Ok(())
}
