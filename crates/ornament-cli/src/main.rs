//! ornament CLI - Modification-aware tRNA scanner
//!
//! Scans genomic sequences for tRNAs and analyzes modification compatibility.

use anyhow::{anyhow, Result};
use clap::{Parser, Subcommand, ValueEnum};
use std::io::Write;
use std::path::Path;

/// Whether the scan output target is a gzip file (`.gz`), which can't be streamed (it needs a
/// final gzip footer written on close).
fn output_is_gz(output: Option<&str>) -> bool {
    output.is_some_and(|p| p.ends_with(".gz"))
}

/// Write `content` to `path`, gzip-compressing transparently when the path ends in `.gz`
/// (matching the scanner's transparent gzip *input* handling).
fn write_output(path: &str, content: &str) -> Result<()> {
    if path.ends_with(".gz") {
        let file = std::fs::File::create(path)?;
        let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        enc.write_all(content.as_bytes())?;
        enc.finish()?;
    } else {
        std::fs::write(path, content)?;
    }
    Ok(())
}

#[derive(Parser)]
#[command(name = "ornament")]
#[command(author, version, about = "Modification-aware tRNA scanner", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

/// Search engine used by the `scan` command.
#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Engine {
    /// Native pure-Rust CM engine (no external binary or C toolchain). Default.
    #[default]
    Native,
    /// External `cmsearch` subprocess (Infernal must be on PATH). Oracle/fallback.
    Cmsearch,
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

        /// Output format (tsv, json, stockholm). `stockholm` writes a per-model
        /// alignment of the hits (native engine only).
        #[arg(short, long, default_value = "tsv")]
        format: String,

        /// Search engine: native (pure Rust, default) or cmsearch (subprocess)
        #[arg(long, value_enum, default_value_t = Engine::Native)]
        engine: Engine,
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

        /// MODOMICS JSON file for modification database (default: built-in)
        #[arg(long)]
        modomics: Option<String>,
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

        /// MODOMICS JSON file for modification database (default: built-in)
        #[arg(long)]
        modomics: Option<String>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Scan {
            input,
            cm,
            output,
            format,
            engine,
        } => {
            use ornament_core::infernal::{
                scan_native_aligned, scan_native_multi, scan_native_multi_streaming,
                write_stockholm, write_tsv, write_tsv_rows, InfernalRunner, TSV_HEADER,
            };

            // E-value reporting threshold (mirrors `cmsearch -E`).
            const E_VALUE: f64 = 1e-5;

            let cm_path = cm.ok_or_else(|| anyhow!("--cm is required"))?;

            // Verify input file exists
            if !Path::new(&input).exists() {
                return Err(anyhow!("Input file not found: {}", input));
            }

            // Verify CM file exists
            if !Path::new(&cm_path).exists() {
                return Err(anyhow!("CM file not found: {}", cm_path));
            }

            // Streaming fast path: native engine + TSV to stdout or a plain file. Hits are written
            // as each model finishes, so a long genome-vs-Rfam run shows progress, survives a
            // crash, and never holds the full hit set in memory. (`.gz` output needs a final gzip
            // footer, so it falls through to the buffered path below, which finalizes correctly.)
            let can_stream =
                engine == Engine::Native && format == "tsv" && !output_is_gz(output.as_deref());
            if can_stream {
                eprintln!(
                    "Scanning {} using {} (native engine, all models, streaming)...",
                    input, cm_path
                );
                let mut w: Box<dyn Write + Send> = match output.as_deref() {
                    Some(p) => Box::new(std::io::BufWriter::new(std::fs::File::create(p)?)),
                    None => Box::new(std::io::BufWriter::new(std::io::stdout())),
                };
                writeln!(w, "{TSV_HEADER}")?;
                let w = std::sync::Mutex::new(w);
                let n = scan_native_multi_streaming(&cm_path, &input, E_VALUE, |hits| {
                    let mut s = String::new();
                    write_tsv_rows(&mut s, hits);
                    let mut g = w.lock().unwrap();
                    let _ = g.write_all(s.as_bytes());
                    let _ = g.flush(); // flush per batch so partial results survive a crash
                })?;
                w.into_inner().unwrap().flush()?;
                eprintln!("Found {n} hits");
                return Ok(());
            }

            let output_str = if format == "stockholm" {
                // Stockholm output keeps each hit's full alignment, so it uses a distinct
                // native scan path. The cmsearch subprocess yields no native alignment.
                if !matches!(engine, Engine::Native) {
                    return Err(anyhow!(
                        "--format stockholm requires --engine native \
                         (the cmsearch subprocess does not produce native alignments)"
                    ));
                }
                eprintln!(
                    "Scanning {} using {} (native engine, Stockholm alignment)...",
                    input, cm_path
                );
                let msas = scan_native_aligned(&cm_path, &input, E_VALUE)?;
                let n_hits: usize = msas.iter().map(|m| m.rows.len()).sum();
                eprintln!("Found {} hits across {} model(s)", n_hits, msas.len());
                write_stockholm(&msas)
            } else {
                let hits = match engine {
                    Engine::Native => {
                        eprintln!(
                            "Scanning {} using {} (native engine, all models)...",
                            input, cm_path
                        );
                        // Runs every model in the .cm file (one, or the whole Rfam collection),
                        // parallelized across the model × record product.
                        scan_native_multi(&cm_path, &input, E_VALUE)?
                    }
                    Engine::Cmsearch => {
                        eprintln!(
                            "Scanning {} for tRNAs using {} (cmsearch subprocess)...",
                            input, cm_path
                        );
                        let runner = InfernalRunner::new()
                            .with_cm(&cm_path)
                            .with_e_value(E_VALUE);
                        runner.cmsearch(&input)?
                    }
                };

                eprintln!("Found {} hits", hits.len());

                match format.as_str() {
                    "json" => serde_json::to_string_pretty(&hits)?,
                    "tsv" => write_tsv(&hits),
                    _ => {
                        return Err(anyhow!(
                            "Unknown format: {}. Use 'tsv', 'json', or 'stockholm'",
                            format
                        ))
                    }
                }
            };

            // Write output
            if let Some(output_path) = output {
                write_output(&output_path, &output_str)?;
                eprintln!("Results written to {}", output_path);
            } else {
                println!("{}", output_str);
            }
        }

        Commands::Analyze {
            input,
            output,
            threshold,
            modomics,
        } => {
            use ornament_core::analysis::{analyze_batch, TRNAHit};

            // Verify input file exists
            if !Path::new(&input).exists() {
                return Err(anyhow!("Input file not found: {}", input));
            }

            eprintln!("Analyzing modification compatibility in {}...", input);
            eprintln!("Threshold: {}", threshold);

            // Read input file (JSON from scan command)
            let content = std::fs::read_to_string(&input)?;
            let hits: Vec<TRNAHit> = serde_json::from_str(&content).map_err(|e| {
                anyhow!(
                    "Failed to parse input JSON: {}. Expected output from 'ornament scan'.",
                    e
                )
            })?;

            eprintln!("Loaded {} tRNA hits", hits.len());

            // Load modification database
            let db = if let Some(modomics_path) = modomics {
                eprintln!("Loading MODOMICS database from {}...", modomics_path);
                ornament_core::modification::ModificationDatabase::from_modomics_file(Path::new(
                    &modomics_path,
                ))
                .map_err(|e| anyhow!("Failed to load MODOMICS file: {}", e))?
            } else {
                ornament_core::modification::ModificationDatabase::eukaryotic()
            };
            let results = analyze_batch(&hits, &db);

            // Filter to odd tRNAs based on threshold
            let odd_results: Vec<_> = results
                .results
                .iter()
                .filter(|r| r.compatibility_score < threshold)
                .collect();

            eprintln!(
                "Found {} odd tRNAs (score < {})",
                odd_results.len(),
                threshold
            );
            eprintln!(
                "Average compatibility: {:.2}%",
                results.average_compatibility * 100.0
            );

            // Format output
            let output_data = serde_json::json!({
                "summary": {
                    "total_trnas": results.total_trnas,
                    "odd_trnas": results.odd_trnas,
                    "average_compatibility": results.average_compatibility,
                    "threshold": threshold
                },
                "odd_trnas": odd_results.iter().map(|r| {
                    serde_json::json!({
                        "id": r.hit.id,
                        "seq_name": r.hit.seq_name,
                        "start": r.hit.start,
                        "end": r.hit.end,
                        "isotype": r.hit.isotype,
                        "anticodon": r.hit.anticodon,
                        "compatibility_score": r.compatibility_score,
                        "incompatibilities": r.incompatibilities.iter().map(|i| {
                            serde_json::json!({
                                "position": i.position.0,
                                "observed_base": i.observed_base.to_char(),
                                "expected_modification": i.expected_mod_name,
                                "severity": format!("{:?}", i.severity)
                            })
                        }).collect::<Vec<_>>()
                    })
                }).collect::<Vec<_>>(),
                "all_results": results.results
            });

            let output_str = serde_json::to_string_pretty(&output_data)?;

            // Write output
            if let Some(output_path) = output {
                write_output(&output_path, &output_str)?;
                eprintln!("Results written to {}", output_path);
            } else {
                println!("{}", output_str);
            }
        }

        Commands::Compare {
            trna,
            modkit,
            output,
        } => {
            use ornament_core::analysis::ModCompatibilityResult;
            use ornament_core::integration::modkit::parse_bedmethyl;

            // Verify input files exist
            if !Path::new(&trna).exists() {
                return Err(anyhow!("tRNA file not found: {}", trna));
            }
            if !Path::new(&modkit).exists() {
                return Err(anyhow!("modkit file not found: {}", modkit));
            }

            eprintln!("Comparing {} with modkit calls from {}...", trna, modkit);

            // Load tRNA analysis results
            let trna_content = std::fs::read_to_string(&trna)?;
            let trna_data: serde_json::Value = serde_json::from_str(&trna_content)?;

            // Extract results from analysis output
            let trna_results: Vec<ModCompatibilityResult> =
                if let Some(results) = trna_data.get("all_results") {
                    serde_json::from_value(results.clone())?
                } else {
                    // Try parsing as direct array of results
                    serde_json::from_str(&trna_content)?
                };

            eprintln!("Loaded {} tRNA results", trna_results.len());

            // Load modkit BedMethyl data
            let modkit_content = std::fs::read_to_string(&modkit)?;
            let modkit_records = parse_bedmethyl(&modkit_content);

            eprintln!("Loaded {} modkit modification calls", modkit_records.len());

            // Build comparison: find modkit calls that overlap with tRNA positions
            let mut comparisons = Vec::new();

            for trna_result in &trna_results {
                let hit = &trna_result.hit;

                // Find modkit records overlapping this tRNA
                let overlapping: Vec<_> = modkit_records
                    .iter()
                    .filter(|r| r.chrom == hit.seq_name && r.start >= hit.start && r.end <= hit.end)
                    .collect();

                if !overlapping.is_empty() {
                    comparisons.push(serde_json::json!({
                        "trna_id": hit.id,
                        "seq_name": hit.seq_name,
                        "start": hit.start,
                        "end": hit.end,
                        "isotype": hit.isotype,
                        "compatibility_score": trna_result.compatibility_score,
                        "is_odd": trna_result.is_odd,
                        "modkit_calls": overlapping.iter().map(|r| {
                            serde_json::json!({
                                "position": r.start,
                                "mod_code": r.mod_code,
                                "strand": r.strand.to_string(),
                                "coverage": r.coverage,
                                "mod_frequency": r.mod_frequency
                            })
                        }).collect::<Vec<_>>(),
                        "expected_incompatibilities": trna_result.incompatibilities.iter().map(|i| {
                            serde_json::json!({
                                "position": i.position.0,
                                "expected_mod": i.expected_mod_name
                            })
                        }).collect::<Vec<_>>()
                    }));
                }
            }

            eprintln!("Found {} tRNAs with modkit overlaps", comparisons.len());

            let output_data = serde_json::json!({
                "summary": {
                    "total_trnas": trna_results.len(),
                    "trnas_with_modkit_calls": comparisons.len(),
                    "total_modkit_records": modkit_records.len()
                },
                "comparisons": comparisons
            });

            let output_str = serde_json::to_string_pretty(&output_data)?;

            if let Some(output_path) = output {
                write_output(&output_path, &output_str)?;
                eprintln!("Results written to {}", output_path);
            } else {
                println!("{}", output_str);
            }
        }

        Commands::Mods {
            position,
            verbose,
            modomics,
        } => {
            let db = if let Some(modomics_path) = modomics {
                eprintln!("Loading MODOMICS database from {}...", modomics_path);
                ornament_core::modification::ModificationDatabase::from_modomics_file(Path::new(
                    &modomics_path,
                ))
                .map_err(|e| anyhow!("Failed to load MODOMICS file: {}", e))?
            } else {
                ornament_core::modification::ModificationDatabase::eukaryotic()
            };

            if let Some(pos) = position {
                let sprinzl = ornament_core::SprinzlPosition(pos.clone());
                let expectations = db.get_expectations(&sprinzl);
                if !expectations.is_empty() {
                    println!("Modifications expected at position {}:", pos);
                    for exp in expectations {
                        for modification in &exp.modifications {
                            if verbose {
                                println!(
                                    "  {} ({}) - {:?} conservation",
                                    modification.name, modification.short_name, exp.conservation
                                );
                            } else {
                                println!("  {}", modification.short_name);
                            }
                        }
                    }
                } else {
                    println!("No modifications expected at position {}", pos);
                }
            } else {
                println!(
                    "Modification database ({} modifications):",
                    db.modifications().len()
                );
                for (name, modification) in db.modifications() {
                    if verbose {
                        println!(
                            "  {} ({}) - derived from {:?}, genomic: {:?}",
                            name,
                            modification.short_name,
                            modification.parent_base,
                            modification.genomic_expectation
                        );
                    } else {
                        println!("  {} ({})", modification.short_name, name);
                    }
                }
            }
        }
    }

    Ok(())
}
