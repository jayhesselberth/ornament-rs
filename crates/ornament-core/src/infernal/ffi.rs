//! Safe FFI wrappers around Infernal, HMMER, and Easel C libraries
//!
//! This module provides safe Rust wrappers around the raw C pointers from infernal-sys.

use anyhow::{anyhow, Result};
use std::ffi::{CStr, CString};
use std::path::Path;
use std::ptr;

use super::CMHit;

// Re-export the raw types for internal use
use infernal_sys::{
    CM_FILE, CM_HIT, CM_PIPELINE, CM_TOPHITS, CM_t, ESL_ALPHABET, ESL_SQ, ESL_SQFILE,
    P7_BG, P7_OPROFILE, P7_SCOREDATA,
};

/// Alphabet type constants from Easel
const ESL_RNA: i32 = 1;

/// CM pipeline mode constants
const CM_SEARCH_SEQS: u32 = 0;

/// Z setby constants
const CM_ZSETBY_SSIINFO: u32 = 0;

/// Safe wrapper around ESL_ALPHABET
pub struct Alphabet {
    ptr: *mut ESL_ALPHABET,
}

impl Alphabet {
    /// Create an RNA alphabet
    pub fn rna() -> Result<Self> {
        unsafe {
            let ptr = infernal_sys::esl_alphabet_Create(ESL_RNA);
            if ptr.is_null() {
                return Err(anyhow!("Failed to create RNA alphabet"));
            }
            Ok(Self { ptr })
        }
    }

    /// Get the raw pointer (for FFI calls)
    pub fn as_ptr(&self) -> *mut ESL_ALPHABET {
        self.ptr
    }
}

impl Drop for Alphabet {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                infernal_sys::esl_alphabet_Destroy(self.ptr);
            }
        }
    }
}

// Alphabet is not Send/Sync due to raw pointer
unsafe impl Send for Alphabet {}

/// Safe wrapper around CM_t (covariance model)
pub struct CovarianceModel {
    ptr: *mut CM_t,
    /// File offset for the CM (needed for cm_Pipeline)
    pub(crate) offset: i64,
}

impl CovarianceModel {
    /// Load a covariance model from a file
    pub fn from_file(path: &Path) -> Result<Self> {
        let path_str = path.to_str().ok_or_else(|| anyhow!("Invalid path"))?;
        let c_path = CString::new(path_str)?;
        let mut errbuf = vec![0u8; 256];

        unsafe {
            let mut cmfp: *mut CM_FILE = ptr::null_mut();

            // Open the CM file
            let status = infernal_sys::cm_file_Open(
                c_path.as_ptr() as *mut i8,
                ptr::null_mut(), // no env
                0,               // don't allow 1.0 format
                &mut cmfp,
                errbuf.as_mut_ptr() as *mut i8,
            );

            if status != 0 || cmfp.is_null() {
                let err_msg = CStr::from_ptr(errbuf.as_ptr() as *const i8)
                    .to_string_lossy()
                    .to_string();
                return Err(anyhow!("Failed to open CM file: {}", err_msg));
            }

            // Read the CM (with fp7 HMM)
            let mut abc: *mut ESL_ALPHABET = ptr::null_mut();
            let mut cm: *mut CM_t = ptr::null_mut();

            let status = infernal_sys::cm_file_Read(
                cmfp,
                1, // read_fp7 = true (read the embedded HMM)
                &mut abc,
                &mut cm,
            );

            // Get the offset before closing
            let offset = if !cmfp.is_null() {
                // The offset is stored in the CM_FILE structure
                // For simplicity, we'll use 0 (single CM file)
                0i64
            } else {
                0i64
            };

            // Close the file
            infernal_sys::cm_file_Close(cmfp);

            if status != 0 || cm.is_null() {
                if !abc.is_null() {
                    infernal_sys::esl_alphabet_Destroy(abc);
                }
                return Err(anyhow!("Failed to read CM from file"));
            }

            // The alphabet is attached to the CM, so we don't need to keep it separately
            // But we should not destroy it here as it's owned by the CM

            Ok(Self { ptr: cm, offset })
        }
    }

    /// Configure the CM for search
    pub fn configure(&mut self) -> Result<()> {
        let mut errbuf = vec![0u8; 256];

        unsafe {
            let status = infernal_sys::cm_Configure(
                self.ptr,
                errbuf.as_mut_ptr() as *mut i8,
                -1, // W_from_cmdline = -1 (use default)
            );

            if status != 0 {
                let err_msg = CStr::from_ptr(errbuf.as_ptr() as *const i8)
                    .to_string_lossy()
                    .to_string();
                return Err(anyhow!("Failed to configure CM: {}", err_msg));
            }
        }

        Ok(())
    }

    /// Get the consensus length of the CM
    pub fn clen(&self) -> i32 {
        unsafe { (*self.ptr).clen }
    }

    /// Get the window length (W) of the CM
    pub fn w(&self) -> i32 {
        unsafe { (*self.ptr).W }
    }

    /// Get the raw pointer
    pub fn as_ptr(&self) -> *mut CM_t {
        self.ptr
    }

    /// Get the alphabet from the CM
    pub fn alphabet(&self) -> *const ESL_ALPHABET {
        unsafe { (*self.ptr).abc }
    }
}

impl Drop for CovarianceModel {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                infernal_sys::FreeCM(self.ptr);
            }
        }
    }
}

unsafe impl Send for CovarianceModel {}

/// Safe wrapper around CM_TOPHITS
pub struct TopHits {
    ptr: *mut CM_TOPHITS,
}

impl TopHits {
    /// Create a new empty top hits container
    pub fn new() -> Result<Self> {
        unsafe {
            let ptr = infernal_sys::cm_tophits_Create();
            if ptr.is_null() {
                return Err(anyhow!("Failed to create TopHits"));
            }
            Ok(Self { ptr })
        }
    }

    /// Get the number of hits
    pub fn len(&self) -> usize {
        unsafe { (*self.ptr).N as usize }
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Sort hits by E-value
    pub fn sort_by_evalue(&mut self) -> Result<()> {
        unsafe {
            let status = infernal_sys::cm_tophits_SortByEvalue(self.ptr);
            if status != 0 {
                return Err(anyhow!("Failed to sort hits by E-value"));
            }
        }
        Ok(())
    }

    /// Convert hits to CMHit structs
    pub fn to_hits(&self, target_name: &str) -> Vec<CMHit> {
        let mut hits = Vec::with_capacity(self.len());

        unsafe {
            let th = &*self.ptr;

            // After sorting, hits are accessible via th.hit array
            // Before sorting, use th.unsrt
            for i in 0..th.N {
                let hit: &CM_HIT = if th.is_sorted_by_evalue != 0 {
                    &*(*th.hit.add(i as usize))
                } else {
                    &*th.unsrt.add(i as usize)
                };

                // Extract hit information
                let name = if !hit.name.is_null() {
                    CStr::from_ptr(hit.name).to_string_lossy().to_string()
                } else {
                    target_name.to_string()
                };

                let (start, end, strand) = if hit.in_rc != 0 {
                    (hit.stop as usize, hit.start as usize, '-')
                } else {
                    (hit.start as usize, hit.stop as usize, '+')
                };

                hits.push(CMHit {
                    target_name: name,
                    target_start: start,
                    target_end: end,
                    strand,
                    query_name: String::new(), // Will be filled in by caller
                    score: hit.score as f64,
                    e_value: hit.evalue,
                    gc_content: 0.0, // Not available directly from CM_HIT
                });
            }
        }

        hits
    }

    /// Get the raw pointer
    pub fn as_ptr(&self) -> *mut CM_TOPHITS {
        self.ptr
    }
}

impl Default for TopHits {
    fn default() -> Self {
        Self::new().expect("Failed to create TopHits")
    }
}

impl Drop for TopHits {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                infernal_sys::cm_tophits_Destroy(self.ptr);
            }
        }
    }
}

unsafe impl Send for TopHits {}

/// Safe wrapper around ESL_SQFILE for reading sequences
pub struct SequenceFile {
    ptr: *mut ESL_SQFILE,
}

impl SequenceFile {
    /// Open a sequence file (FASTA format)
    pub fn open(path: &Path, alphabet: &Alphabet) -> Result<Self> {
        let path_str = path.to_str().ok_or_else(|| anyhow!("Invalid path"))?;
        let c_path = CString::new(path_str)?;

        unsafe {
            let mut sqfp: *mut ESL_SQFILE = ptr::null_mut();

            // Open with format autodetection (format = 0)
            let status = infernal_sys::esl_sqfile_Open(
                c_path.as_ptr() as *mut i8,
                0, // format = unknown (autodetect)
                ptr::null_mut(), // env
                &mut sqfp,
            );

            if status != 0 || sqfp.is_null() {
                return Err(anyhow!("Failed to open sequence file: {}", path_str));
            }

            // Set digital mode with the alphabet
            let status = infernal_sys::esl_sqfile_SetDigital(sqfp, alphabet.as_ptr());
            if status != 0 {
                infernal_sys::esl_sqfile_Close(sqfp);
                return Err(anyhow!("Failed to set digital mode for sequence file"));
            }

            Ok(Self { ptr: sqfp })
        }
    }

    /// Get the raw pointer
    pub fn as_ptr(&self) -> *mut ESL_SQFILE {
        self.ptr
    }
}

impl Drop for SequenceFile {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                infernal_sys::esl_sqfile_Close(self.ptr);
            }
        }
    }
}

unsafe impl Send for SequenceFile {}

/// Safe wrapper around ESL_SQ (sequence)
pub struct Sequence {
    ptr: *mut ESL_SQ,
}

impl Sequence {
    /// Create a digital sequence for the given alphabet
    pub fn create_digital(alphabet: &Alphabet) -> Result<Self> {
        unsafe {
            let ptr = infernal_sys::esl_sq_CreateDigital(alphabet.as_ptr());
            if ptr.is_null() {
                return Err(anyhow!("Failed to create digital sequence"));
            }
            Ok(Self { ptr })
        }
    }

    /// Get the sequence name
    pub fn name(&self) -> String {
        unsafe {
            let name_ptr = (*self.ptr).name;
            if name_ptr.is_null() {
                String::new()
            } else {
                CStr::from_ptr(name_ptr).to_string_lossy().to_string()
            }
        }
    }

    /// Get the sequence length
    pub fn len(&self) -> i64 {
        unsafe { (*self.ptr).n }
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Reuse the sequence buffer for reading next sequence
    pub fn reuse(&mut self) -> Result<()> {
        unsafe {
            let status = infernal_sys::esl_sq_Reuse(self.ptr);
            if status != 0 {
                return Err(anyhow!("Failed to reuse sequence buffer"));
            }
        }
        Ok(())
    }

    /// Get the raw pointer
    pub fn as_ptr(&self) -> *mut ESL_SQ {
        self.ptr
    }
}

impl Drop for Sequence {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe {
                infernal_sys::esl_sq_Destroy(self.ptr);
            }
        }
    }
}

unsafe impl Send for Sequence {}

/// Read the next sequence from a file into the sequence buffer
pub fn read_sequence(sqfp: &SequenceFile, sq: &mut Sequence) -> Result<bool> {
    unsafe {
        let status = infernal_sys::esl_sqio_Read(sqfp.as_ptr(), sq.as_ptr());

        // eslEOF = 1 means end of file
        if status == 1 {
            return Ok(false);
        }

        if status != 0 {
            return Err(anyhow!("Error reading sequence (status: {})", status));
        }

        Ok(true)
    }
}

/// HMM filter components needed for cm_Pipeline
/// These are extracted from the CM's embedded HMM
pub struct HmmFilter {
    pub(crate) om: *mut P7_OPROFILE,
    pub(crate) bg: *mut P7_BG,
    pub(crate) msvdata: *mut P7_SCOREDATA,
    /// E-value parameters (array of floats)
    pub(crate) p7_evparam: Vec<f32>,
}

impl HmmFilter {
    /// Create HMM filter from a configured CM
    ///
    /// The CM must have been configured and have an embedded HMM (fp7).
    pub fn from_cm(cm: &CovarianceModel, alphabet: &Alphabet) -> Result<Self> {
        unsafe {
            let cm_ptr = cm.as_ptr();
            let cm_ref = &*cm_ptr;

            // Check if CM has embedded HMM
            if cm_ref.fp7.is_null() {
                return Err(anyhow!("CM does not have embedded HMM (fp7)"));
            }

            // Get HMM model length from the embedded HMM
            let hmm_m = (*cm_ref.fp7).M;

            // Create background model
            let bg = infernal_sys::p7_bg_Create(alphabet.as_ptr());
            if bg.is_null() {
                return Err(anyhow!("Failed to create P7 background model"));
            }

            // Create optimized profile from the CM's fp7 HMM
            let om = infernal_sys::p7_oprofile_Create(hmm_m, alphabet.as_ptr());
            if om.is_null() {
                infernal_sys::p7_bg_Destroy(bg);
                return Err(anyhow!("Failed to create P7 optimized profile"));
            }

            // Create a standard profile (needed for configuration and conversion)
            let gm = infernal_sys::p7_profile_Create(hmm_m, alphabet.as_ptr());
            if gm.is_null() {
                infernal_sys::p7_oprofile_Destroy(om);
                infernal_sys::p7_bg_Destroy(bg);
                return Err(anyhow!("Failed to create P7 profile"));
            }

            // Configure the profile (L=100 is a typical target length, mode=2 is p7_LOCAL)
            let status = infernal_sys::p7_ProfileConfig(cm_ref.fp7, bg, gm, 100, 2);
            if status != 0 {
                infernal_sys::p7_profile_Destroy(gm);
                infernal_sys::p7_oprofile_Destroy(om);
                infernal_sys::p7_bg_Destroy(bg);
                return Err(anyhow!("Failed to configure P7 profile"));
            }

            // Convert to optimized profile
            let status = infernal_sys::p7_oprofile_Convert(gm, om);
            if status != 0 {
                infernal_sys::p7_profile_Destroy(gm);
                infernal_sys::p7_oprofile_Destroy(om);
                infernal_sys::p7_bg_Destroy(bg);
                return Err(anyhow!("Failed to convert to optimized profile"));
            }

            // Create MSV data (needs both om and gm)
            let msvdata = infernal_sys::p7_hmm_ScoreDataCreate(om, gm);
            if msvdata.is_null() {
                infernal_sys::p7_profile_Destroy(gm);
                infernal_sys::p7_oprofile_Destroy(om);
                infernal_sys::p7_bg_Destroy(bg);
                return Err(anyhow!("Failed to create MSV score data"));
            }

            // Copy E-value parameters from CM's fp7_evparam
            // The array has 8 elements in CM, but we use 6 for p7
            let mut p7_evparam = vec![0.0f32; 8];
            for i in 0..8 {
                p7_evparam[i] = cm_ref.fp7_evparam[i];
            }

            // Clean up the standard profile (we only need the optimized one)
            infernal_sys::p7_profile_Destroy(gm);

            Ok(Self {
                om,
                bg,
                msvdata,
                p7_evparam,
            })
        }
    }
}

impl Drop for HmmFilter {
    fn drop(&mut self) {
        unsafe {
            if !self.msvdata.is_null() {
                infernal_sys::p7_hmm_ScoreDataDestroy(self.msvdata);
            }
            if !self.om.is_null() {
                infernal_sys::p7_oprofile_Destroy(self.om);
            }
            if !self.bg.is_null() {
                infernal_sys::p7_bg_Destroy(self.bg);
            }
        }
    }
}

unsafe impl Send for HmmFilter {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alphabet_creation() {
        let abc = Alphabet::rna();
        assert!(abc.is_ok());
    }

    #[test]
    fn test_tophits_creation() {
        let th = TopHits::new();
        assert!(th.is_ok());
        let th = th.unwrap();
        assert_eq!(th.len(), 0);
        assert!(th.is_empty());
    }
}
