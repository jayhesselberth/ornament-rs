//! Glocal CYK alignment with traceback (ported from the CYK align + `cm_alignT` in
//! `src/cm_dpalign.c`).
//!
//! Given a hit window `i0..j0`, this aligns the whole subsequence to the model (ROOT->END),
//! recovers the optimal parse, and—via the [`EmitMap`]—reports which **consensus column**
//! each residue occupies. That per-residue consensus mapping is exactly what the Sprinzl
//! mapper needs, and is where modification annotation attaches.

use ornament_alphabet::Dsq;

use crate::config::IMPOSSIBLE;
use crate::emit::Oesc;
use crate::emitmap::EmitMap;
use crate::model::{emits_right, n_emit, st, Cm};

/// One aligned residue: its 1-based sequence position and, if it aligns to a match column,
/// the consensus column (1..clen). Insert residues have `consensus = None`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignedResidue {
    pub seq_pos: usize,
    pub consensus: Option<usize>,
    /// True for a paired (MP) emission, useful for reconstructing SS_cons.
    pub paired: bool,
}

/// Result of aligning a hit window.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// CYK alignment score (bits).
    pub score: f32,
    /// Residues in sequence order, each tagged with its consensus column (if a match).
    pub residues: Vec<AlignedResidue>,
}

impl Alignment {
    /// Consensus column for each sequence position, as a `Vec` aligned to `i0..j0`
    /// (`None` for inserts). Convenience for the Sprinzl mapper.
    pub fn consensus_columns(&self) -> Vec<Option<usize>> {
        self.residues.iter().map(|r| r.consensus).collect()
    }
}

#[inline]
fn add(a: f32, b: f32) -> f32 {
    if a <= IMPOSSIBLE || b <= IMPOSSIBLE {
        IMPOSSIBLE
    } else {
        a + b
    }
}

/// Align the hit window `dsq[i0..=j0]` (1-based, inclusive) to the model in glocal mode.
/// Returns the CYK score and the per-residue consensus mapping.
pub fn align_glocal(cm: &Cm, dsq: &[Dsq], i0: usize, j0: usize, emap: &EmitMap) -> Alignment {
    align_glocal_impl(cm, dsq, i0, j0, emap, None)
}

/// The CYK-optimal glocal parse as a list of `(v, i, j, d)` for every visited CM state (`i`, `j`
/// **window-local**, i.e. `1..=L`; `d = j - i + 1`). Used to validate that HMM bands contain the
/// true parse (the no-drop guarantee).
#[cfg(test)]
pub(crate) fn align_glocal_parse(
    cm: &Cm,
    dsq: &[Dsq],
    i0: usize,
    j0: usize,
    emap: &EmitMap,
) -> Vec<(usize, usize, usize, usize)> {
    let mut parse = Vec::new();
    align_glocal_impl(cm, dsq, i0, j0, emap, Some(&mut parse));
    parse
}

fn align_glocal_impl(
    cm: &Cm,
    dsq: &[Dsq],
    i0: usize,
    j0: usize,
    emap: &EmitMap,
    mut parse: Option<&mut Vec<(usize, usize, usize, usize)>>,
) -> Alignment {
    assert!(
        i0 >= 1 && j0 >= i0 && j0 + 1 < dsq.len(),
        "window out of range"
    );
    let l = j0 - i0 + 1;
    let m = cm.m;
    let oesc = Oesc::build(cm);
    let kp = oesc.kp;

    // Full DP matrix alpha[v][j][d], j,d in 0..=L. (Hit windows are short.)
    let sj = l + 1; // stride over d
    let sv = sj * sj; // stride over j*d
    let idx = |v: usize, j: usize, d: usize| v * sv + j * sj + d;
    let mut alpha = vec![IMPOSSIBLE; m * sv];
    let mut yshadow = vec![-1i32; m * sv]; // chosen child yoffset (non-B)
    let mut kshadow = vec![-1i32; m * sv]; // chosen bifurcation split (B)

    // Local sequence access: local position p in 1..=L -> code.
    let res = |p: usize| dsq[i0 + p - 1] as usize;

    for j in 0..=l {
        for v in (0..m).rev() {
            let sttype = cm.sttype[v];
            if sttype == st::E {
                alpha[idx(v, j, 0)] = 0.0;
                continue;
            }
            let sd = n_emit(sttype);
            let right = emits_right(sttype);

            if sttype == st::B {
                let wch = cm.cfirst[v] as usize;
                let ych = cm.cnum[v] as usize;
                for d in 0..=j {
                    let mut best = IMPOSSIBLE;
                    let mut bk = -1i32;
                    for k in 0..=d {
                        let cand = add(alpha[idx(wch, j - k, d - k)], alpha[idx(ych, j, k)]);
                        if cand > best {
                            best = cand;
                            bk = k as i32;
                        }
                    }
                    alpha[idx(v, j, d)] = best;
                    kshadow[idx(v, j, d)] = bk;
                }
                continue;
            }

            let ych = cm.cfirst[v] as usize;
            let cnum = cm.cnum[v] as usize;
            let tsc = &cm.tsc[v];
            let single = oesc.single[v].as_deref();
            let pairtab = oesc.pair[v].as_deref();

            for d in 0..=j {
                if d < sd {
                    continue; // emitters need >= sd residues
                }
                let jp = if right { j - 1 } else { j };
                let cd = d - sd;
                let mut best = IMPOSSIBLE;
                let mut by = -1i32;
                for yoff in 0..cnum {
                    let cand = add(alpha[idx(ych + yoff, jp, cd)], tsc[yoff]);
                    if cand > best {
                        best = cand;
                        by = yoff as i32;
                    }
                }
                let val = match sttype {
                    st::ML | st::IL => add(best, single.unwrap()[res(j - d + 1)]),
                    st::MR | st::IR => add(best, single.unwrap()[res(j)]),
                    st::MP => add(best, pairtab.unwrap()[res(j - d + 1) * kp + res(j)]),
                    _ => best, // D, S
                };
                alpha[idx(v, j, d)] = val;
                yshadow[idx(v, j, d)] = by;
            }
        }
    }

    let score = alpha[idx(0, l, l)];

    // --- Traceback from (ROOT_S, j=L, d=L). ---
    let mut residues: Vec<AlignedResidue> = Vec::new();
    let mut stack: Vec<(usize, usize, usize)> = vec![(0, l, l)];
    while let Some((v, j, d)) = stack.pop() {
        let sttype = cm.sttype[v];
        let nd = cm.ndidx[v];
        if let Some(p) = parse.as_deref_mut() {
            // window-local i = j - d + 1 (i > j when d == 0, e.g. E states).
            p.push((v, j + 1 - d, j, d));
        }
        match sttype {
            st::ML | st::IL => {
                let p = j - d + 1; // local left position
                residues.push(AlignedResidue {
                    seq_pos: i0 + p - 1,
                    consensus: if sttype == st::ML {
                        Some(emap.lpos[nd] as usize)
                    } else {
                        None
                    },
                    paired: false,
                });
            }
            st::MR | st::IR => {
                residues.push(AlignedResidue {
                    seq_pos: i0 + j - 1,
                    consensus: if sttype == st::MR {
                        Some(emap.rpos[nd] as usize)
                    } else {
                        None
                    },
                    paired: false,
                });
            }
            st::MP => {
                let p = j - d + 1;
                residues.push(AlignedResidue {
                    seq_pos: i0 + p - 1,
                    consensus: Some(emap.lpos[nd] as usize),
                    paired: true,
                });
                residues.push(AlignedResidue {
                    seq_pos: i0 + j - 1,
                    consensus: Some(emap.rpos[nd] as usize),
                    paired: true,
                });
            }
            _ => {}
        }

        if sttype == st::E {
            continue;
        }
        if sttype == st::B {
            let k = kshadow[idx(v, j, d)];
            if k < 0 {
                continue;
            }
            let k = k as usize;
            let wch = cm.cfirst[v] as usize;
            let ych = cm.cnum[v] as usize;
            stack.push((wch, j - k, d - k)); // left fragment (BEGL subtree)
            stack.push((ych, j, k)); // right fragment (BEGR subtree)
        } else {
            let by = yshadow[idx(v, j, d)];
            if by < 0 {
                continue;
            }
            let child = cm.cfirst[v] as usize + by as usize;
            let sd = n_emit(sttype);
            let jp = if emits_right(sttype) { j - 1 } else { j };
            stack.push((child, jp, d - sd));
        }
    }

    residues.sort_by_key(|r| r.seq_pos);
    Alignment { score, residues }
}
