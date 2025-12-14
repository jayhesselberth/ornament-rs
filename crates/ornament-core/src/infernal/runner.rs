//! Infernal command runner
//!
//! Executes cmsearch as subprocess for tRNA detection.

use std::path::Path;
use std::process::Command;
use anyhow::{anyhow, Result};

use super::parser::parse_tblout;
use super::CMHit;

/// Runner for Infernal commands
pub struct InfernalRunner {
    cm_path: Option<String>,
    /// E-value threshold for reporting hits
    e_value: f64,
    /// Number of CPUs to use
    cpu: usize,
}

impl InfernalRunner {
    pub fn new() -> Self {
        Self {
            cm_path: None,
            e_value: 1e-5,
            cpu: num_cpus(),
        }
    }

    /// Set the covariance model file path
    pub fn with_cm<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.cm_path = Some(path.as_ref().to_string_lossy().to_string());
        self
    }

    /// Set E-value threshold
    pub fn with_e_value(mut self, e: f64) -> Self {
        self.e_value = e;
        self
    }

    /// Set number of CPUs
    pub fn with_cpu(mut self, cpu: usize) -> Self {
        self.cpu = cpu;
        self
    }

    /// Run cmsearch on a FASTA file using subprocess
    pub fn cmsearch<P: AsRef<Path>>(&self, fasta: P) -> Result<Vec<CMHit>> {
        let cm_path = self
            .cm_path
            .as_ref()
            .ok_or_else(|| anyhow!("No covariance model specified"))?;

        let fasta_path = fasta.as_ref();
        if !fasta_path.exists() {
            return Err(anyhow!("FASTA file not found: {}", fasta_path.display()));
        }

        // Run cmsearch with tabular output to stdout
        let output = Command::new("cmsearch")
            .arg("--tblout")
            .arg("/dev/stdout")  // Write tabular output to stdout
            .arg("-o")
            .arg("/dev/null")    // Suppress alignment output
            .arg("--cpu")
            .arg(self.cpu.to_string())
            .arg("-E")
            .arg(self.e_value.to_string())
            .arg(cm_path)
            .arg(fasta_path)
            .output()
            .map_err(|e| {
                if e.kind() == std::io::ErrorKind::NotFound {
                    anyhow!("cmsearch not found. Please install Infernal: http://eddylab.org/infernal/")
                } else {
                    anyhow!("Failed to run cmsearch: {}", e)
                }
            })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(anyhow!("cmsearch failed: {}", stderr));
        }

        // Parse the tabular output from stdout
        let tblout_content = String::from_utf8_lossy(&output.stdout);
        let hits = parse_tblout(&tblout_content);

        Ok(hits)
    }
}

impl Default for InfernalRunner {
    fn default() -> Self {
        Self::new()
    }
}

/// Get number of available CPUs
fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_runner_builder() {
        let runner = InfernalRunner::new()
            .with_cm("test.cm")
            .with_e_value(1e-10)
            .with_cpu(4);

        assert_eq!(runner.cm_path, Some("test.cm".to_string()));
        assert_eq!(runner.e_value, 1e-10);
        assert_eq!(runner.cpu, 4);
    }
}
