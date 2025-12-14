//! Infernal command runner
//!
//! Executes cmsearch/cmscan as subprocesses for tRNA detection.

use std::path::Path;
use anyhow::Result;

/// Runner for Infernal commands
pub struct InfernalRunner {
    cm_path: Option<String>,
}

impl InfernalRunner {
    pub fn new() -> Self {
        Self { cm_path: None }
    }

    pub fn with_cm<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.cm_path = Some(path.as_ref().to_string_lossy().to_string());
        self
    }

    /// Run cmsearch on a FASTA file
    pub fn cmsearch<P: AsRef<Path>>(&self, _fasta: P) -> Result<Vec<super::CMHit>> {
        // TODO: Implement subprocess call to cmsearch
        Ok(vec![])
    }
}

impl Default for InfernalRunner {
    fn default() -> Self {
        Self::new()
    }
}
