//! Output helpers shared across commands.

use anyhow::Result;
use std::io::Write;

/// Whether the output target is a gzip file (`.gz`), which can't be streamed (it needs a
/// final gzip footer written on close).
pub fn output_is_gz(output: Option<&str>) -> bool {
    output.is_some_and(|p| p.ends_with(".gz"))
}

/// Write `content` to `path`, gzip-compressing transparently when the path ends in `.gz`
/// (matching the scanner's transparent gzip *input* handling).
pub fn write_output(path: &str, content: &str) -> Result<()> {
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

/// Write a fully-rendered `content` string to `path` (gzip-aware) or stdout when `path` is
/// `None`. Logs the destination for file output.
pub fn emit(output: Option<&str>, content: &str) -> Result<()> {
    match output {
        Some(p) => {
            write_output(p, content)?;
            tracing::info!("Results written to {}", crate::style::path(p));
        }
        None => {
            println!("{content}");
        }
    }
    Ok(())
}
