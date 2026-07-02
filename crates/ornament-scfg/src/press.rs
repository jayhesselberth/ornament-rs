//! Prebuilt CM database ("press"), analogous to Infernal's `cmpress`.
//!
//! Every `cm_pipeline_search` call rebuilds a model's p7 filter cascade and QDB bands; a many-model
//! scan additionally reparses the ASCII `.cm(.gz)` and reconfigures scores on every run. [`press_cm_file`]
//! does that per-model preparation once and serializes it to a binary `<db>.orm` sidecar;
//! [`load_pressed`] reads it back, skipping the parse + configure + O(#B·W²) band computation.
//!
//! What the file stores: the local-configured scan model, the glocal alignment model, both models'
//! QDB bands, and the parsed p7 filter HMM — the expensive, *portable* inputs. The SIMD-striped
//! filter profiles are **not** stored (their byte layout is target-CPU-dependent); they are cheaply
//! rebuilt on load via [`PreparedFilters::from_parts`]. The alphabet is stored once as a compact
//! [`AlphabetDescriptor`] and reattached to every model on load.
//!
//! The codec is [`postcard`] (a compact, non-self-describing serde format), so any change to the
//! serialized structs must bump [`PRESS_VERSION`]; [`load_pressed`] hard-rejects a magic/version
//! mismatch rather than decoding garbage.

use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

use ornament_alphabet::{Alphabet, AlphabetDescriptor};
use ornament_hmm::P7Hmm;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cmfile::parse_cm_records_file;
use crate::config::{configure_local, configure_scores};
use crate::model::Cm;
use crate::pipeline::{parse_filter_hmm, PreparedFilters};
use crate::qdb::{calc_qdb_bands, QdbBands};
use crate::InfernalError;

/// Magic bytes identifying a pressed CM database.
pub const PRESS_MAGIC: [u8; 8] = *b"ORMPRESS";
/// On-disk format version. Bump on any change to the serialized structs; [`load_pressed`] rejects
/// files whose version differs (the format is not self-describing, so a mismatch can't be decoded safely).
pub const PRESS_VERSION: u32 = 1;
/// Sidecar file extension appended to a `.cm` path when pressing (`Rfam.cm` → `Rfam.cm.orm`).
pub const PRESS_EXT: &str = "orm";

/// File-level header: identifies the format and records how the bands were built.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PressHeader {
    pub magic: [u8; 8],
    pub version: u32,
    /// `ornament-scfg` crate version that wrote the file (diagnostic; a mismatch warns, not fails).
    pub crate_version: String,
    /// QDB tail-loss `beta` the stored bands were computed with.
    pub beta: f64,
    pub n_models: u32,
}

/// One model's fully-prepared, record-independent search state, as persisted.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PressedModel {
    /// Local-configured model (the `cmsearch`/`cmscan` default) used for scanning.
    pub cm: Cm,
    /// Glocal-configured copy used only for per-hit model-span alignment.
    pub align_cm: Cm,
    /// QDB bands for the scan DP, computed on the local `cm`.
    pub scan_bands: QdbBands,
    /// QDB bands for the alignment DP, computed on the glocal `align_cm`.
    pub align_bands: QdbBands,
    /// Parsed p7 filter HMM (from the CM's embedded fp7 text).
    pub fp7: P7Hmm,
    /// Max hit width (`cm.w`), the scan window ceiling.
    pub w_max: usize,
}

impl PressedModel {
    /// Rebuild this model's [`PreparedFilters`] (striped p7 profiles + the stored scan bands),
    /// skipping the fp7 parse and the QDB band computation those inputs already cover.
    pub fn prepared_filters(&self) -> Result<PreparedFilters, InfernalError> {
        PreparedFilters::from_parts(&self.cm, self.fp7.clone(), self.scan_bands.clone())
    }
}

/// The whole serialized collection. Private: callers go through [`press_cm_file`]/[`load_pressed`].
#[derive(Debug, Serialize, Deserialize)]
struct PressedCollection {
    header: PressHeader,
    /// The alphabet, stored once and reattached to every model on load.
    alphabet: AlphabetDescriptor,
    models: Vec<PressedModel>,
}

/// A loaded pressed database: the shared alphabet plus every prepared model (with alphabets reattached).
#[derive(Debug)]
pub struct PressedDb {
    pub abc: Arc<Alphabet>,
    pub models: Vec<PressedModel>,
}

/// The default sidecar path for a `.cm` file: append `.orm` to the full name (`Rfam.cm` → `Rfam.cm.orm`).
pub fn default_press_path(cm_path: &Path) -> std::path::PathBuf {
    let mut s = cm_path.as_os_str().to_os_string();
    s.push(".");
    s.push(PRESS_EXT);
    std::path::PathBuf::from(s)
}

/// Prepare every model in `cm_path` and write the pressed DB to `out_path`. Returns the model count.
///
/// Prep mirrors the per-model setup in the scan pipeline: a glocal copy for alignment, local config
/// for scanning, both models' QDB bands, and the parsed filter HMM — run in parallel across models.
/// Refuses to overwrite an existing `out_path` unless `force`. If `out_path` ends in `.gz` the output
/// is gzip-compressed.
pub fn press_cm_file(cm_path: &Path, out_path: &Path, force: bool) -> Result<usize, InfernalError> {
    if out_path.exists() && !force {
        return Err(InfernalError::Io(format!(
            "refusing to overwrite existing {} (use --force)",
            out_path.display()
        )));
    }

    let models = parse_cm_records_file(cm_path)?;
    if models.is_empty() {
        return Err(InfernalError::Parse(format!(
            "no CM models found in {}",
            cm_path.display()
        )));
    }
    // Every model parsed from one `.cm` file shares the same alphabet; capture its recipe once.
    let alphabet = models[0].abc.descriptor();
    let beta = QdbBands::DEFAULT_BETA;

    let pressed: Vec<PressedModel> = models
        .into_par_iter()
        .map(|mut cm| -> Result<PressedModel, InfernalError> {
            let mut align_cm = cm.clone();
            configure_scores(&mut align_cm);
            configure_local(&mut cm); // cmsearch/cmscan default (local mode)
            let scan_bands = calc_qdb_bands(&cm, beta);
            let align_bands = calc_qdb_bands(&align_cm, beta);
            let fp7 = parse_filter_hmm(&cm)?;
            let w_max = cm.w as usize;
            Ok(PressedModel {
                cm,
                align_cm,
                scan_bands,
                align_bands,
                fp7,
                w_max,
            })
        })
        .collect::<Result<Vec<_>, InfernalError>>()?;

    let n = pressed.len();
    let collection = PressedCollection {
        header: PressHeader {
            magic: PRESS_MAGIC,
            version: PRESS_VERSION,
            crate_version: env!("CARGO_PKG_VERSION").to_string(),
            beta,
            n_models: n as u32,
        },
        alphabet,
        models: pressed,
    };

    let bytes = postcard::to_allocvec(&collection)
        .map_err(|e| InfernalError::Io(format!("failed to encode pressed DB: {e}")))?;
    write_bytes(out_path, &bytes)?;
    Ok(n)
}

/// Load a pressed DB written by [`press_cm_file`], reattaching the shared alphabet to every model.
///
/// Errors on a magic or [`PRESS_VERSION`] mismatch (the file was written by an incompatible build).
/// A differing `crate_version` is tolerated (the format version already gates compatibility).
pub fn load_pressed(path: &Path) -> Result<PressedDb, InfernalError> {
    let bytes = read_bytes(path)?;
    let mut collection: PressedCollection = postcard::from_bytes(&bytes).map_err(|e| {
        InfernalError::Parse(format!(
            "{} is not a valid pressed DB (or was written by an incompatible version): {e}",
            path.display()
        ))
    })?;

    if collection.header.magic != PRESS_MAGIC {
        return Err(InfernalError::Parse(format!(
            "{} is not a pressed CM database (bad magic)",
            path.display()
        )));
    }
    if collection.header.version != PRESS_VERSION {
        return Err(InfernalError::Parse(format!(
            "{} has pressed-format version {}, but this build expects {} — re-run `ornament press`",
            path.display(),
            collection.header.version,
            PRESS_VERSION
        )));
    }

    // Rebuild the one shared alphabet and reattach it to every (scan + align) model, replacing the
    // placeholder each `Cm` deserialized with. Sharing one `Arc` keeps memory flat across the DB.
    let abc = Arc::new(collection.alphabet.build()?);
    for m in &mut collection.models {
        m.cm.abc = Arc::clone(&abc);
        m.align_cm.abc = Arc::clone(&abc);
    }

    Ok(PressedDb {
        abc,
        models: collection.models,
    })
}

/// Read a peek of `path`'s header without decoding the whole DB — returns `None` if it isn't a
/// pressed file (bad magic / undecodable), so callers can auto-detect a stale/foreign sidecar.
pub fn peek_header(path: &Path) -> Option<PressHeader> {
    let bytes = read_bytes(path).ok()?;
    let collection: PressedCollection = postcard::from_bytes(&bytes).ok()?;
    if collection.header.magic == PRESS_MAGIC {
        Some(collection.header)
    } else {
        None
    }
}

/// Write `bytes` to `path`, gzip-compressing when the path ends in `.gz` (mirrors the CLI's output I/O).
fn write_bytes(path: &Path, bytes: &[u8]) -> Result<(), InfernalError> {
    let file = std::fs::File::create(path)
        .map_err(|e| InfernalError::Io(format!("failed to create {}: {e}", path.display())))?;
    if is_gz(path) {
        let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        enc.write_all(bytes)
            .and_then(|_| enc.finish().map(|_| ()))
            .map_err(|e| InfernalError::Io(format!("failed to write {}: {e}", path.display())))
    } else {
        let mut w = std::io::BufWriter::new(file);
        w.write_all(bytes)
            .and_then(|_| w.flush())
            .map_err(|e| InfernalError::Io(format!("failed to write {}: {e}", path.display())))
    }
}

/// Read all bytes from `path`, transparently gunzipping a gzip stream (magic-sniffed, extension-independent).
fn read_bytes(path: &Path) -> Result<Vec<u8>, InfernalError> {
    let raw = std::fs::read(path)
        .map_err(|e| InfernalError::Io(format!("failed to read {}: {e}", path.display())))?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .map_err(|e| {
                InfernalError::Io(format!("failed to decompress {}: {e}", path.display()))
            })?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

fn is_gz(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()) == Some("gz")
}
