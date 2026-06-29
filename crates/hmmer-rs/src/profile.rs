//! Search profile (`P7_PROFILE`) built from a [`P7Hmm`] — a faithful port of
//! `p7_ProfileConfig` + `p7_ReconfigLength` (HMMER `modelconfig.c`) for **multihit local**
//! mode (`p7_LOCAL`), the configuration Infernal uses for its CM filter.
//!
//! All scores are natural-log log-odds (nats), matching the generic Forward/Viterbi/MSV DP.

use easel_rs::Alphabet;

use crate::hmm::{tr, P7Hmm};

/// Transition-score indices into `tsc[k*8 + .]` (HMMER `p7P_*`).
pub mod p7p {
    pub const MM: usize = 0;
    pub const IM: usize = 1;
    pub const DM: usize = 2;
    pub const BM: usize = 3;
    pub const MD: usize = 4;
    pub const DD: usize = 5;
    pub const MI: usize = 6;
    pub const II: usize = 7;
    pub const NTRANS: usize = 8;
}

/// Special-state indices for `xsc[state][LOOP|MOVE]` (HMMER `p7P_E/N/J/C/B`).
pub mod xstate {
    pub const E: usize = 0;
    pub const N: usize = 1;
    pub const J: usize = 2;
    pub const C: usize = 3;
    pub const B: usize = 4;
    pub const LOOP: usize = 0;
    pub const MOVE: usize = 1;
}

/// A configured Plan7 search profile (multihit local).
#[derive(Debug, Clone)]
pub struct P7Profile {
    pub m: usize,
    pub kp: usize,
    /// Transition scores `tsc[k*8 + type]`, k in 0..=M.
    pub tsc: Vec<f32>,
    /// Match emission scores `msc[x][k]`, x in 0..Kp, k in 0..=M.
    pub msc: Vec<Vec<f32>>,
    /// Insert emission scores `isc[x][k]`.
    pub isc: Vec<Vec<f32>>,
    /// Special transitions `xsc[state][LOOP|MOVE]`.
    pub xsc: [[f32; 2]; 5],
    pub nj: f32,
    pub l: usize,
}

fn ln(x: f32) -> f32 {
    x.ln() // ln(0) = -inf, matching HMMER's log(0)
}

/// Expected match-state occupancy (`p7_hmm_CalculateOccupancy`, mocc only).
fn occupancy(hmm: &P7Hmm) -> Vec<f32> {
    let m = hmm.m;
    let mut occ = vec![0.0f32; m + 1];
    if m >= 1 {
        occ[1] = hmm.t[0][tr::MI] + hmm.t[0][tr::MM];
    }
    for k in 2..=m {
        occ[k] = occ[k - 1] * (hmm.t[k - 1][tr::MM] + hmm.t[k - 1][tr::MI])
            + (1.0 - occ[k - 1]) * hmm.t[k - 1][tr::DM];
    }
    occ
}

impl P7Profile {
    /// Configure a multihit-local profile for target length `l` (`p7_ProfileConfig` with
    /// `mode = p7_LOCAL`, then `p7_ReconfigLength`). `abc` supplies the alphabet for the
    /// degenerate emission expansion; `bg` is the background frequency vector (uniform 1/K
    /// for RNA, as `p7_bg_Create`).
    pub fn config_local(hmm: &P7Hmm, abc: &Alphabet, bg: &[f32], l: usize) -> P7Profile {
        let m = hmm.m;
        let k = hmm.k;
        let kp = abc.kp;
        let ninf = f32::NEG_INFINITY;

        let mut tsc = vec![ninf; (m + 1) * p7p::NTRANS];

        // Local entry (BM): occ[k] / sum_i occ[i]*(M-i+1); stored at index k-1.
        let occ = occupancy(hmm);
        let z: f32 = (1..=m).map(|i| occ[i] * (m - i + 1) as f32).sum();
        for kk in 1..=m {
            tsc[(kk - 1) * p7p::NTRANS + p7p::BM] = ln(occ[kk] / z);
        }

        // Node transitions (k = 1..M-1).
        for kk in 1..m {
            let t = &hmm.t[kk];
            let base = kk * p7p::NTRANS;
            tsc[base + p7p::MM] = ln(t[tr::MM]);
            tsc[base + p7p::MI] = ln(t[tr::MI]);
            tsc[base + p7p::MD] = ln(t[tr::MD]);
            tsc[base + p7p::IM] = ln(t[tr::IM]);
            tsc[base + p7p::II] = ln(t[tr::II]);
            tsc[base + p7p::DM] = ln(t[tr::DM]);
            tsc[base + p7p::DD] = ln(t[tr::DD]);
        }

        // Match emission scores. Canonical: log(mat[k][x]/bg[x]); degenerate via the
        // alphabet's marginalization (favg_score), which reduces to the canonical value.
        let mut msc = vec![vec![ninf; m + 1]; kp];
        let mut isc = vec![vec![ninf; m + 1]; kp];
        for kk in 1..=m {
            let mut canon = vec![ninf; k];
            for (x, c) in canon.iter_mut().enumerate() {
                *c = ln(hmm.mat[kk][x] / bg[x]);
            }
            for (x, row) in msc.iter_mut().enumerate() {
                row[kk] = if abc.is_residue(x as u8) {
                    abc.favg_score(x as u8, &canon)
                } else {
                    ninf
                };
            }
        }
        // Insert emissions hardwired to 0 (background) for k < M; I_M impossible.
        for (x, row) in isc.iter_mut().enumerate() {
            if abc.is_residue(x as u8) {
                row[1..m].fill(0.0);
            }
        }

        // E state: multihit local -> MOVE = LOOP = -log(2); nj = 1.
        let mut xsc = [[ninf; 2]; 5];
        let log2 = std::f32::consts::LN_2;
        xsc[xstate::E][xstate::MOVE] = -log2;
        xsc[xstate::E][xstate::LOOP] = -log2;
        let nj = 1.0f32;

        let mut prof = P7Profile {
            m,
            kp,
            tsc,
            msc,
            isc,
            xsc,
            nj,
            l: 0,
        };
        prof.reconfig_length(l);
        prof
    }

    /// `p7_ReconfigLength`: set the N/C/J length-model transitions for target length `l`.
    pub fn reconfig_length(&mut self, l: usize) {
        let pmove = (2.0 + self.nj) / (l as f32 + 2.0 + self.nj);
        let ploop = 1.0 - pmove;
        for st in [xstate::N, xstate::C, xstate::J] {
            self.xsc[st][xstate::LOOP] = ploop.ln();
            self.xsc[st][xstate::MOVE] = pmove.ln();
        }
        self.l = l;
    }

    #[inline]
    pub fn tsc(&self, ttype: usize, k: usize) -> f32 {
        self.tsc[k * p7p::NTRANS + ttype]
    }
}

/// Default background frequencies (`p7_bg_Create`): uniform `1/K` for nucleic alphabets.
pub fn bg_freqs(k: usize) -> Vec<f32> {
    vec![1.0 / k as f32; k]
}

/// Null-model score for a length-`L` sequence (`p7_bg_NullOne` with `p7_bg_SetLength`):
/// `L*log(p1) + log(1-p1)`, `p1 = L/(L+1)`. In nats.
pub fn null_one(l: usize) -> f32 {
    let p1 = l as f32 / (l as f32 + 1.0);
    l as f32 * p1.ln() + (1.0 - p1).ln()
}
