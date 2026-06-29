//! Optimized emission scores (`oesc`) indexed by digital code, shared by the scanning and
//! alignment DPs.
//!
//! Singlet emitters get a length-`Kp` table; MP states get `Kp*Kp`. Degenerate codes are
//! marginalized and modified-base symbols resolve to their parent base — this is the single
//! point where the modification-aware alphabet flows into the DP.

use easel_rs::Dsq;

use crate::config::IMPOSSIBLE;
use crate::model::{st, Cm};

pub(crate) struct Oesc {
    pub kp: usize,
    pub single: Vec<Option<Vec<f32>>>, // singlet emitters: length Kp
    pub pair: Vec<Option<Vec<f32>>>,   // MP: length Kp*Kp
}

impl Oesc {
    pub(crate) fn build(cm: &Cm) -> Oesc {
        let abc = &cm.abc;
        let kp = abc.kp;
        let k = abc.k;
        let mut single = vec![None; cm.m];
        let mut pair = vec![None; cm.m];

        for v in 0..cm.m {
            match cm.sttype[v] {
                st::ML | st::MR | st::IL | st::IR => {
                    let esc = &cm.esc[v];
                    let row: Vec<f32> = (0..kp)
                        .map(|a| {
                            if abc.is_residue(a as Dsq) {
                                abc.favg_score(a as Dsq, esc)
                            } else {
                                IMPOSSIBLE
                            }
                        })
                        .collect();
                    single[v] = Some(row);
                }
                st::MP => {
                    let esc = &cm.esc[v]; // length k*k
                    let mut row = vec![IMPOSSIBLE; kp * kp];
                    for a in 0..kp {
                        if !abc.is_residue(a as Dsq) {
                            continue;
                        }
                        for b in 0..kp {
                            if !abc.is_residue(b as Dsq) {
                                continue;
                            }
                            let mut sum = 0.0f32;
                            let mut n = 0u32;
                            for x in 0..k {
                                if !abc.degen[a][x] {
                                    continue;
                                }
                                for y in 0..k {
                                    if abc.degen[b][y] {
                                        sum += esc[x * k + y];
                                        n += 1;
                                    }
                                }
                            }
                            row[a * kp + b] = if n > 0 { sum / n as f32 } else { IMPOSSIBLE };
                        }
                    }
                    pair[v] = Some(row);
                }
                _ => {}
            }
        }
        Oesc { kp, single, pair }
    }
}
