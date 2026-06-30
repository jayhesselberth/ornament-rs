//! `ornament trna` — the tRNA-modification analysis application, built on the general engine.

use std::path::Path;

use anyhow::{anyhow, Result};
use clap::{Args, Subcommand};
use tabled::builder::Builder;
use tabled::settings::Style;

use ornament_core::modification::ModificationDatabase;

use super::require_file;
use crate::io::emit;
use crate::style;

#[derive(Subcommand, Debug)]
pub enum TrnaCommand {
    /// Analyze modification compatibility of tRNA sequences
    Analyze(AnalyzeArgs),
    /// Compare tRNA analysis with modkit modification calls
    Compare(CompareArgs),
    /// Show the modification database
    Mods(ModsArgs),
}

pub fn run(cmd: TrnaCommand) -> Result<()> {
    match cmd {
        TrnaCommand::Analyze(a) => analyze(a),
        TrnaCommand::Compare(a) => compare(a),
        TrnaCommand::Mods(a) => mods(a),
    }
}

/// Load the modification database from a MODOMICS file, or the built-in eukaryotic default.
fn load_db(modomics: Option<&str>) -> Result<ModificationDatabase> {
    if let Some(p) = modomics {
        tracing::info!("Loading MODOMICS database from {}", style::path(p));
        ModificationDatabase::from_modomics_file(Path::new(p))
            .map_err(|e| anyhow!("Failed to load MODOMICS file: {e}"))
    } else {
        Ok(ModificationDatabase::eukaryotic())
    }
}

// ----------------------------------------------------------------------------- analyze

#[derive(Args, Debug)]
pub struct AnalyzeArgs {
    /// Input tRNA hits (JSON from `ornament scan`/`search`)
    #[arg(short, long)]
    pub input: String,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<String>,

    /// Compatibility score threshold for "odd" tRNAs
    #[arg(short, long, default_value = "0.8")]
    pub threshold: f64,

    /// MODOMICS JSON file for modification database (default: built-in)
    #[arg(long)]
    pub modomics: Option<String>,
}

fn analyze(args: AnalyzeArgs) -> Result<()> {
    use ornament_core::analysis::{analyze_batch, TRNAHit};

    require_file(&args.input, "Input file")?;

    tracing::info!(
        "{} modification compatibility in {} ({} {})",
        style::action("Analyzing"),
        style::path(&args.input),
        style::label("threshold:"),
        style::value(args.threshold),
    );

    let content = std::fs::read_to_string(&args.input)?;
    let hits: Vec<TRNAHit> = serde_json::from_str(&content).map_err(|e| {
        anyhow!("Failed to parse input JSON: {e}. Expected output from 'ornament scan'.")
    })?;
    tracing::info!("Loaded {} tRNA hits", style::count(hits.len()));

    let db = load_db(args.modomics.as_deref())?;
    let results = analyze_batch(&hits, &db);

    let odd_results: Vec<_> = results
        .results
        .iter()
        .filter(|r| r.compatibility_score < args.threshold)
        .collect();

    tracing::info!(
        "Found {} odd tRNAs (score < {})",
        style::count(odd_results.len()),
        args.threshold
    );
    tracing::info!(
        "Average compatibility: {}",
        style::percentage(format!("{:.2}%", results.average_compatibility * 100.0))
    );

    let output_data = serde_json::json!({
        "summary": {
            "total_trnas": results.total_trnas,
            "odd_trnas": results.odd_trnas,
            "average_compatibility": results.average_compatibility,
            "threshold": args.threshold
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

    emit(
        args.output.as_deref(),
        &serde_json::to_string_pretty(&output_data)?,
    )
}

// ----------------------------------------------------------------------------- compare

#[derive(Args, Debug)]
pub struct CompareArgs {
    /// tRNA analysis results (JSON from `ornament trna analyze`)
    #[arg(short, long)]
    pub trna: String,

    /// modkit BedMethyl file
    #[arg(short, long)]
    pub modkit: String,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<String>,
}

fn compare(args: CompareArgs) -> Result<()> {
    use ornament_core::analysis::ModCompatibilityResult;
    use ornament_core::integration::modkit::parse_bedmethyl;

    require_file(&args.trna, "tRNA file")?;
    require_file(&args.modkit, "modkit file")?;

    tracing::info!(
        "{} {} with modkit calls from {}",
        style::action("Comparing"),
        style::path(&args.trna),
        style::path(&args.modkit),
    );

    let trna_content = std::fs::read_to_string(&args.trna)?;
    let trna_data: serde_json::Value = serde_json::from_str(&trna_content)?;

    let trna_results: Vec<ModCompatibilityResult> =
        if let Some(results) = trna_data.get("all_results") {
            serde_json::from_value(results.clone())?
        } else {
            serde_json::from_str(&trna_content)?
        };
    tracing::info!("Loaded {} tRNA results", style::count(trna_results.len()));

    let modkit_content = std::fs::read_to_string(&args.modkit)?;
    let modkit_records = parse_bedmethyl(&modkit_content);
    tracing::info!(
        "Loaded {} modkit modification calls",
        style::count(modkit_records.len())
    );

    let mut comparisons = Vec::new();
    for trna_result in &trna_results {
        let hit = &trna_result.hit;
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

    tracing::info!(
        "Found {} tRNAs with modkit overlaps",
        style::count(comparisons.len())
    );

    let output_data = serde_json::json!({
        "summary": {
            "total_trnas": trna_results.len(),
            "trnas_with_modkit_calls": comparisons.len(),
            "total_modkit_records": modkit_records.len()
        },
        "comparisons": comparisons
    });

    emit(
        args.output.as_deref(),
        &serde_json::to_string_pretty(&output_data)?,
    )
}

// ----------------------------------------------------------------------------- mods

#[derive(Args, Debug)]
pub struct ModsArgs {
    /// Filter by Sprinzl position (e.g., "34", "55")
    #[arg(short, long)]
    pub position: Option<String>,

    /// Show detailed information (named `--long` to avoid clashing with the global `-v`)
    #[arg(short = 'l', long)]
    pub long: bool,

    /// MODOMICS JSON file for modification database (default: built-in)
    #[arg(long)]
    pub modomics: Option<String>,
}

fn mods(args: ModsArgs) -> Result<()> {
    let db = load_db(args.modomics.as_deref())?;

    if let Some(pos) = args.position {
        let sprinzl = ornament_core::SprinzlPosition(pos.clone());
        let expectations = db.get_expectations(&sprinzl);
        if expectations.is_empty() {
            println!("No modifications expected at position {pos}");
            return Ok(());
        }
        eprintln!(
            "{}",
            style::header(format!("Modifications expected at position {pos}"))
        );
        let mut builder = Builder::default();
        if args.long {
            builder.push_record(["modification", "short", "conservation"]);
            for exp in expectations {
                for m in &exp.modifications {
                    builder.push_record([
                        m.name.clone(),
                        m.short_name.clone(),
                        format!("{:?}", exp.conservation),
                    ]);
                }
            }
        } else {
            builder.push_record(["short", "modification"]);
            for exp in expectations {
                for m in &exp.modifications {
                    builder.push_record([m.short_name.clone(), m.name.clone()]);
                }
            }
        }
        println!("{}", builder.build().with(Style::rounded()));
        return Ok(());
    }

    let modifications = db.modifications();
    eprintln!(
        "{} ({} modifications)",
        style::header("Modification database"),
        modifications.len()
    );
    let mut builder = Builder::default();
    if args.long {
        builder.push_record(["short", "name", "parent", "genomic"]);
        for (name, m) in modifications {
            builder.push_record([
                m.short_name.clone(),
                name.clone(),
                format!("{:?}", m.parent_base),
                format!("{:?}", m.genomic_expectation),
            ]);
        }
    } else {
        builder.push_record(["short", "name"]);
        for (name, m) in modifications {
            builder.push_record([m.short_name.clone(), name.clone()]);
        }
    }
    println!("{}", builder.build().with(Style::rounded()));
    Ok(())
}
