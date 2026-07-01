//! Model statistics — the native `cmstat` analog.
//!
//! Reads a `.cm` (one model or a whole collection) and extracts the metadata
//! `cmstat` reports: name/accession, consensus length, node/state counts, the
//! training-set size, Pfam-style GA/TC/NC thresholds, and whether the model
//! carries calibrated E-value parameters. The heavy lifting is just field access
//! on [`ornament_scfg::Cm`]; this module turns it into a serializable record the
//! CLI renders (as a table) without depending on `ornament-scfg` directly.

use anyhow::{anyhow, Result};
use std::path::Path;

use ornament_scfg::model::exp_mode;
use ornament_scfg::{parse_cm_records_file, Cm};

/// One model's summary statistics (a row of `ornament stat` output).
#[derive(Debug, Clone, serde::Serialize)]
pub struct CmStat {
    /// 1-based index within the `.cm` file (matches `cmstat`'s `idx` column).
    pub idx: usize,
    pub name: String,
    pub accession: Option<String>,
    pub description: Option<String>,
    /// Number of states (`Cm::m`).
    pub states: usize,
    pub nodes: usize,
    /// Consensus length (`2·MATP + MATL + MATR`).
    pub clen: usize,
    /// Maximum hit length the model was built for.
    pub w: i32,
    /// Number of training sequences.
    pub nseq: i32,
    /// Effective sequence count after entropy weighting.
    pub eff_nseq: f32,
    pub ga: Option<f32>,
    pub tc: Option<f32>,
    pub nc: Option<f32>,
    /// True when the model carries calibrated E-value parameters (glocal-CYK tail valid).
    pub calibrated: bool,
    /// Glocal-CYK exponential tail `(lambda, mu)` when calibrated.
    pub glocal_cyk: Option<(f64, f64)>,
}

impl CmStat {
    fn from_cm(idx: usize, cm: &Cm) -> Self {
        let gc = &cm.exp[exp_mode::CM_GC];
        let calibrated = gc.is_valid;
        CmStat {
            idx: idx + 1,
            name: cm.name.clone(),
            accession: cm.acc.clone(),
            description: cm.desc.clone(),
            states: cm.m,
            nodes: cm.nodes,
            clen: cm.clen,
            w: cm.w,
            nseq: cm.nseq,
            eff_nseq: cm.eff_nseq,
            ga: cm.ga,
            tc: cm.tc,
            nc: cm.nc,
            calibrated,
            glocal_cyk: calibrated.then_some((gc.lambda, gc.mu_extrap)),
        }
    }
}

/// Parse every model in `cm_path` (gzip-transparent) and return its summary statistics,
/// in file order. The native equivalent of running `cmstat` on a `.cm` file.
pub fn cmstat<P: AsRef<Path>>(cm_path: P) -> Result<Vec<CmStat>> {
    let cm_path = cm_path.as_ref();
    let models = parse_cm_records_file(cm_path)
        .map_err(|e| anyhow!("failed to parse CM file {}: {e}", cm_path.display()))?;
    Ok(models
        .iter()
        .enumerate()
        .map(|(i, cm)| CmStat::from_cm(i, cm))
        .collect())
}
