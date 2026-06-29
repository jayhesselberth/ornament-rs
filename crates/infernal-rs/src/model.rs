//! The covariance model (`CM_t`) data structure, ported faithfully from Infernal.
//!
//! Field names and semantics mirror `src/infernal.h`. Emission/transition arrays are
//! sized from the alphabet's `K` (4 for RNA) rather than literal 4/16, so the model is
//! alphabet-size-generic and ready for `K>4` extended-emission discovery models (Phase 8)
//! without touching the DP core.

use std::sync::Arc;

use easel_rs::Alphabet;

/// Maximum number of state-to-state connections per state (`MAXCONNECT`).
pub const MAXCONNECT: usize = 6;

/// State types (`*_st`), values matching Infernal exactly.
pub mod st {
    pub const D: u8 = 0; // delete
    pub const MP: u8 = 1; // match-pair
    pub const ML: u8 = 2; // match-left
    pub const MR: u8 = 3; // match-right
    pub const IL: u8 = 4; // insert-left
    pub const IR: u8 = 5; // insert-right
    pub const S: u8 = 6; // start
    pub const E: u8 = 7; // end
    pub const B: u8 = 8; // bifurcation
    pub const EL: u8 = 9; // local end
}

/// Node types (`*_nd`), values matching Infernal exactly.
pub mod nd {
    pub const BIF: u8 = 0;
    pub const MATP: u8 = 1;
    pub const MATL: u8 = 2;
    pub const MATR: u8 = 3;
    pub const BEGL: u8 = 4;
    pub const BEGR: u8 = 5;
    pub const ROOT: u8 = 6;
    pub const END: u8 = 7;
}

/// Unique state identifiers (`stid`), values matching Infernal exactly.
pub mod stid {
    pub const ROOT_S: u8 = 0;
    pub const ROOT_IL: u8 = 1;
    pub const ROOT_IR: u8 = 2;
    pub const BEGL_S: u8 = 3;
    pub const BEGR_S: u8 = 4;
    pub const BEGR_IL: u8 = 5;
    pub const MATP_MP: u8 = 6;
    pub const MATP_ML: u8 = 7;
    pub const MATP_MR: u8 = 8;
    pub const MATP_D: u8 = 9;
    pub const MATP_IL: u8 = 10;
    pub const MATP_IR: u8 = 11;
    pub const MATL_ML: u8 = 12;
    pub const MATL_D: u8 = 13;
    pub const MATL_IL: u8 = 14;
    pub const MATR_MR: u8 = 15;
    pub const MATR_D: u8 = 16;
    pub const MATR_IR: u8 = 17;
    pub const END_E: u8 = 18;
    pub const BIF_B: u8 = 19;
    pub const END_EL: u8 = 20;
}

/// E-value exponential-tail search modes (`EXP_*`), values matching Infernal exactly.
pub mod exp_mode {
    pub const CM_GC: usize = 0; // glocal CYK
    pub const CM_GI: usize = 1; // glocal Inside
    pub const CM_LC: usize = 2; // local CYK
    pub const CM_LI: usize = 3; // local Inside
    pub const NMODES: usize = 4;
}

/// `StateCode`: state-type token -> code.
pub fn state_code(s: &str) -> Option<u8> {
    Some(match s {
        "D" => st::D,
        "MP" => st::MP,
        "ML" => st::ML,
        "MR" => st::MR,
        "IL" => st::IL,
        "IR" => st::IR,
        "S" => st::S,
        "E" => st::E,
        "B" => st::B,
        "EL" => st::EL,
        _ => return None,
    })
}

/// `NodeCode`: node-type token -> code.
pub fn node_code(s: &str) -> Option<u8> {
    Some(match s {
        "BIF" => nd::BIF,
        "MATP" => nd::MATP,
        "MATL" => nd::MATL,
        "MATR" => nd::MATR,
        "BEGL" => nd::BEGL,
        "BEGR" => nd::BEGR,
        "ROOT" => nd::ROOT,
        "END" => nd::END,
        _ => return None,
    })
}

/// `DeriveUniqueStateCode`: (node type, state type) -> unique state id.
pub fn derive_unique_state_code(ndtype: u8, sttype: u8) -> Option<u8> {
    use stid::*;
    Some(match (ndtype, sttype) {
        (nd::BIF, st::B) => BIF_B,
        (nd::MATP, st::D) => MATP_D,
        (nd::MATP, st::MP) => MATP_MP,
        (nd::MATP, st::ML) => MATP_ML,
        (nd::MATP, st::MR) => MATP_MR,
        (nd::MATP, st::IL) => MATP_IL,
        (nd::MATP, st::IR) => MATP_IR,
        (nd::MATL, st::D) => MATL_D,
        (nd::MATL, st::ML) => MATL_ML,
        (nd::MATL, st::IL) => MATL_IL,
        (nd::MATR, st::D) => MATR_D,
        (nd::MATR, st::MR) => MATR_MR,
        (nd::MATR, st::IR) => MATR_IR,
        (nd::BEGL, st::S) => BEGL_S,
        (nd::BEGR, st::S) => BEGR_S,
        (nd::BEGR, st::IL) => BEGR_IL,
        (nd::ROOT, st::S) => ROOT_S,
        (nd::ROOT, st::IL) => ROOT_IL,
        (nd::ROOT, st::IR) => ROOT_IR,
        (nd::END, st::E) => END_E,
        _ => return None,
    })
}

/// Number of residues a state emits: 2 for MP, 1 for ML/MR/IL/IR, 0 otherwise
/// (`StateLeftDelta`/`StateRightDelta` semantics summed).
pub fn n_emit(sttype: u8) -> usize {
    match sttype {
        st::MP => 2,
        st::ML | st::MR | st::IL | st::IR => 1,
        _ => 0,
    }
}

/// Does this state type emit on the left (`StateLeftDelta == 1`)?
pub fn emits_left(sttype: u8) -> bool {
    matches!(sttype, st::MP | st::ML | st::IL)
}

/// Does this state type emit on the right (`StateRightDelta == 1`)?
pub fn emits_right(sttype: u8) -> bool {
    matches!(sttype, st::MP | st::MR | st::IR)
}

/// Exponential-tail fit describing the random-score distribution for one search mode
/// (`ExpInfo_t`). Read from the `.cm` `ECM*` lines; used for CM E-values.
#[derive(Debug, Clone, Copy)]
pub struct ExpInfo {
    pub lambda: f64,
    pub mu_extrap: f64,
    pub mu_orig: f64,
    pub dbsize: f64,
    pub nrandhits: i64,
    pub tailp: f64,
    pub is_valid: bool,
}

impl Default for ExpInfo {
    fn default() -> Self {
        ExpInfo {
            lambda: 0.0,
            mu_extrap: 0.0,
            mu_orig: 0.0,
            dbsize: 0.0,
            nrandhits: 0,
            tailp: 0.0,
            is_valid: false,
        }
    }
}

/// A covariance model. Probabilities (`t`, `e`) are the canonical stored form, parsed via
/// `2^score · null`; log-odds scores (`tsc`, `esc`) are computed by `config`.
#[derive(Debug, Clone)]
pub struct Cm {
    // Annotation / metadata.
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    pub rf: Option<String>,        // 1..clen reference line, if CMH_RF
    pub consensus: Option<String>, // 1..clen consensus residues, if CMH_CONS
    pub map: Option<Vec<i32>>,     // 1..clen alignment-column map, if CMH_MAP (map[0]=0)

    pub nseq: i32,
    pub eff_nseq: f32,

    /// Null (background) residue probabilities, length `abc.K`.
    pub null: Vec<f32>,

    // Core topology.
    pub m: usize,     // number of states
    pub clen: usize,  // consensus length (2*MATP + MATL + MATR)
    pub nodes: usize, // number of nodes
    pub w: i32,       // max hit length

    pub sttype: Vec<u8>, // [M]
    pub ndidx: Vec<usize>, // [M] node each state belongs to
    pub stid: Vec<u8>,   // [M] unique state id
    pub cfirst: Vec<i32>, // [M] first child (or BEGL child for B_st)
    pub cnum: Vec<i32>,  // [M] # connections (or BEGR child index for B_st)
    pub plast: Vec<i32>, // [M] first parent
    pub pnum: Vec<i32>,  // [M] # parents

    pub nodemap: Vec<usize>, // [nodes] first state of each node
    pub ndtype: Vec<u8>,     // [nodes]
    pub node_annot: Vec<NodeAnnot>, // [nodes] raw map/cons/rf columns

    /// QDB bands (read but not yet used by the non-banded DP).
    pub dmin1: Vec<i32>,
    pub dmin2: Vec<i32>,
    pub dmax1: Vec<i32>,
    pub dmax2: Vec<i32>,

    /// Transition probabilities, `[M][0..cnum]`.
    pub t: Vec<Vec<f32>>,
    /// Emission probabilities, `[M]`: length `K` for singlet emitters, `K*K` for MP.
    pub e: Vec<Vec<f32>>,

    /// Log-odds scores, filled by `config`.
    pub tsc: Vec<Vec<f32>>,
    pub esc: Vec<Vec<f32>>,
    /// Local begin/end probabilities, filled by `config`.
    pub begin: Vec<f32>,
    pub end: Vec<f32>,

    pub pbegin: f32,
    pub pend: f32,
    pub el_selfsc: f32,
    pub ga: Option<f32>,
    pub tc: Option<f32>,
    pub nc: Option<f32>,

    /// E-value exponential tails, indexed by [`exp_mode`].
    pub exp: [ExpInfo; exp_mode::NMODES],
    /// p7 glocal-forward E-value params (mu, lambda) from `EFP7GF`.
    pub fp7_gfmu: Option<f64>,
    pub fp7_gflambda: Option<f64>,
    /// Raw text of the embedded HMMER3 filter HMM, retained for Phase 5.
    pub fp7_text: Option<String>,

    pub flags: CmFlags,
    pub abc: Arc<Alphabet>,
}

/// Header status flags (subset of Infernal's `CMH_*`).
#[derive(Debug, Clone, Copy, Default)]
pub struct CmFlags {
    pub rf: bool,
    pub cons: bool,
    pub map: bool,
    pub ga: bool,
    pub tc: bool,
    pub nc: bool,
}

/// Per-node annotation columns from the `.cm` node lines (alignment-column map,
/// consensus residues, and reference characters for the left/right positions). Retained
/// so the consensus / emit-map can be assembled when alignment lands (Phase 4).
#[derive(Debug, Clone, Default)]
pub struct NodeAnnot {
    pub lmap: Option<i32>,
    pub rmap: Option<i32>,
    pub lcons: Option<char>,
    pub rcons: Option<char>,
    pub lrf: Option<char>,
    pub rrf: Option<char>,
}

impl Cm {
    /// True if state `v` is a bifurcation (`B_st`).
    #[inline]
    pub fn is_bif(&self, v: usize) -> bool {
        self.sttype[v] == st::B
    }

    /// True if state `v` emits any residue.
    #[inline]
    pub fn is_emitter(&self, v: usize) -> bool {
        n_emit(self.sttype[v]) > 0
    }

    /// Alphabet size `K`.
    #[inline]
    pub fn k(&self) -> usize {
        self.abc.k
    }
}
