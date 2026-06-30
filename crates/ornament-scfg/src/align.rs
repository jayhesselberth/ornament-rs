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
use crate::qdb::QdbBands;

/// One aligned residue: its 1-based sequence position and, if it aligns to a match column,
/// the consensus column (1..clen). Insert residues have `consensus = None`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignedResidue {
    pub seq_pos: usize,
    pub consensus: Option<usize>,
    /// True for a paired (MP) emission, useful for reconstructing SS_cons.
    pub paired: bool,
    /// For an insert residue (`consensus == None`), the consensus column this insert
    /// follows in the multiple alignment: IL after `lpos[nd]`, IR after `rpos[nd]-1`
    /// (Infernal's `Parsetrees2Alignment` insert-column convention). `None` for matches.
    pub insert_after: Option<usize>,
    /// True for an IR (right) insert, which is right-flushed *after* the IL inserts that
    /// share the same gap (`insert_after`). `false` for IL inserts and match residues.
    pub insert_right: bool,
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

/// QDB-banded glocal alignment: same result as [`align_glocal`] but the DP restricts each state's
/// subsequence length `d` to its query-dependent band `[dmin[v], dmax[v]]`, so the matrix is
/// `O(m·window·band)` instead of the dense `O(m·window²)`. Since the QDB bands contain the optimal
/// parse (they're "safe"), the recovered alignment is **identical** to the unbanded one — this just
/// bounds memory so aligning a long hit on a large model (e.g. rRNA) doesn't blow up. `bands` must
/// be the model's [`QdbBands`].
pub fn align_glocal_banded(
    cm: &Cm,
    dsq: &[Dsq],
    i0: usize,
    j0: usize,
    emap: &EmitMap,
    bands: &QdbBands,
) -> Alignment {
    assert!(
        i0 >= 1 && j0 >= i0 && j0 + 1 < dsq.len(),
        "window out of range"
    );
    let l = j0 - i0 + 1;
    let m = cm.m;
    let oesc = Oesc::build(cm);
    let kp = oesc.kp;

    // Per-state d-band `[dn, dx]` (QDB, clamped to `[0, l]`); each state's block is `(l+1) × bw`
    // cells, laid out at packed prefix offsets `base[v]`.
    let mut dn = vec![0usize; m];
    let mut dx = vec![0usize; m];
    let mut bw = vec![0usize; m];
    let mut base = vec![0usize; m + 1];
    for v in 0..m {
        let lo = bands.dmin[v].min(l);
        let hi = bands.dmax[v].min(l);
        dn[v] = lo;
        dx[v] = hi;
        bw[v] = if hi >= lo { hi - lo + 1 } else { 0 };
        base[v + 1] = base[v] + (l + 1) * bw[v];
    }
    let total = base[m];
    let mut alpha = vec![IMPOSSIBLE; total];
    let mut yshadow = vec![-1i32; total];
    let mut kshadow = vec![-1i32; total];

    let res = |p: usize| dsq[i0 + p - 1] as usize;
    // Banded cell index, or `None` if `d` is out of state `v`'s band (`j` is always `0..=l`).
    let idx = |v: usize, j: usize, d: usize| -> Option<usize> {
        if bw[v] == 0 || d < dn[v] || d > dx[v] {
            None
        } else {
            Some(base[v] + j * bw[v] + (d - dn[v]))
        }
    };

    for j in 0..=l {
        for v in (0..m).rev() {
            if bw[v] == 0 {
                continue; // unreachable state in the bands
            }
            let sttype = cm.sttype[v];
            if sttype == st::E {
                if let Some(ix) = idx(v, j, 0) {
                    alpha[ix] = 0.0; // E: d == 0 only
                }
                continue;
            }
            let sd = n_emit(sttype);
            let right = emits_right(sttype);

            if sttype == st::B {
                let wch = cm.cfirst[v] as usize;
                let ych = cm.cnum[v] as usize;
                for d in dn[v]..=dx[v].min(j) {
                    let mut best = IMPOSSIBLE;
                    let mut bk = -1i32;
                    for k in 0..=d {
                        let left = idx(wch, j - k, d - k).map_or(IMPOSSIBLE, |ix| alpha[ix]);
                        let rgt = idx(ych, j, k).map_or(IMPOSSIBLE, |ix| alpha[ix]);
                        let cand = add(left, rgt);
                        if cand > best {
                            best = cand;
                            bk = k as i32;
                        }
                    }
                    let ix = idx(v, j, d).unwrap();
                    alpha[ix] = best;
                    kshadow[ix] = bk;
                }
                continue;
            }

            let ych = cm.cfirst[v] as usize;
            let cnum = cm.cnum[v] as usize;
            let tsc = &cm.tsc[v];
            let single = oesc.single[v].as_deref();
            let pairtab = oesc.pair[v].as_deref();

            for d in dn[v].max(sd)..=dx[v].min(j) {
                let jp = if right { j - 1 } else { j };
                let cd = d - sd;
                let mut best = IMPOSSIBLE;
                let mut by = -1i32;
                #[allow(clippy::needless_range_loop)] // yoff indexes both the child base and tsc
                for yoff in 0..cnum {
                    let cand = add(
                        idx(ych + yoff, jp, cd).map_or(IMPOSSIBLE, |ix| alpha[ix]),
                        tsc[yoff],
                    );
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
                let ix = idx(v, j, d).unwrap();
                alpha[ix] = val;
                yshadow[ix] = by;
            }
        }
    }

    let score = idx(0, l, l).map_or(IMPOSSIBLE, |ix| alpha[ix]);

    // --- Traceback from (ROOT_S, j=L, d=L), identical to the dense version but banded-indexed. ---
    let mut residues: Vec<AlignedResidue> = Vec::new();
    let mut stack: Vec<(usize, usize, usize)> = vec![(0, l, l)];
    while let Some((v, j, d)) = stack.pop() {
        let sttype = cm.sttype[v];
        let nd = cm.ndidx[v];
        emit_residue(cm, emap, &mut residues, v, sttype, i0, j, d);

        if sttype == st::E {
            continue;
        }
        let Some(ix) = idx(v, j, d) else { continue };
        if sttype == st::B {
            let k = kshadow[ix];
            if k < 0 {
                continue;
            }
            let k = k as usize;
            stack.push((cm.cfirst[v] as usize, j - k, d - k)); // left (BEGL)
            stack.push((cm.cnum[v] as usize, j, k)); // right (BEGR)
        } else {
            let by = yshadow[ix];
            if by < 0 {
                continue;
            }
            let child = cm.cfirst[v] as usize + by as usize;
            let jp = if emits_right(sttype) { j - 1 } else { j };
            stack.push((child, jp, d - n_emit(sttype)));
        }
        let _ = nd;
    }

    residues.sort_by_key(|r| r.seq_pos);
    Alignment { score, residues }
}

/// Push the aligned residue(s) emitted by state `v` at `(j, d)` (shared by the dense and banded
/// tracebacks). `i0` maps window-local positions back to sequence coordinates.
#[allow(clippy::too_many_arguments)]
fn emit_residue(
    cm: &Cm,
    emap: &EmitMap,
    residues: &mut Vec<AlignedResidue>,
    v: usize,
    sttype: u8,
    i0: usize,
    j: usize,
    d: usize,
) {
    let nd = cm.ndidx[v];
    match sttype {
        st::ML | st::IL => {
            let p = j - d + 1;
            let is_match = sttype == st::ML;
            residues.push(AlignedResidue {
                seq_pos: i0 + p - 1,
                consensus: is_match.then(|| emap.lpos[nd] as usize),
                paired: false,
                insert_after: (!is_match).then(|| emap.lpos[nd] as usize),
                insert_right: false,
            });
        }
        st::MR | st::IR => {
            let is_match = sttype == st::MR;
            residues.push(AlignedResidue {
                seq_pos: i0 + j - 1,
                consensus: is_match.then(|| emap.rpos[nd] as usize),
                paired: false,
                insert_after: (!is_match).then(|| (emap.rpos[nd] - 1) as usize),
                insert_right: !is_match,
            });
        }
        st::MP => {
            let p = j - d + 1;
            residues.push(AlignedResidue {
                seq_pos: i0 + p - 1,
                consensus: Some(emap.lpos[nd] as usize),
                paired: true,
                insert_after: None,
                insert_right: false,
            });
            residues.push(AlignedResidue {
                seq_pos: i0 + j - 1,
                consensus: Some(emap.rpos[nd] as usize),
                paired: true,
                insert_after: None,
                insert_right: false,
            });
        }
        _ => {}
    }
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
                let is_match = sttype == st::ML;
                residues.push(AlignedResidue {
                    seq_pos: i0 + p - 1,
                    consensus: is_match.then(|| emap.lpos[nd] as usize),
                    paired: false,
                    // IL inserts follow the node's left consensus column.
                    insert_after: (!is_match).then(|| emap.lpos[nd] as usize),
                    insert_right: false,
                });
            }
            st::MR | st::IR => {
                let is_match = sttype == st::MR;
                residues.push(AlignedResidue {
                    seq_pos: i0 + j - 1,
                    consensus: is_match.then(|| emap.rpos[nd] as usize),
                    paired: false,
                    // IR inserts precede the node's right column, i.e. follow `rpos-1`,
                    // and are right-flushed within the gap (Parsetrees2Alignment).
                    insert_after: (!is_match).then(|| (emap.rpos[nd] - 1) as usize),
                    insert_right: !is_match,
                });
            }
            st::MP => {
                let p = j - d + 1;
                residues.push(AlignedResidue {
                    seq_pos: i0 + p - 1,
                    consensus: Some(emap.lpos[nd] as usize),
                    paired: true,
                    insert_after: None,
                    insert_right: false,
                });
                residues.push(AlignedResidue {
                    seq_pos: i0 + j - 1,
                    consensus: Some(emap.rpos[nd] as usize),
                    paired: true,
                    insert_after: None,
                    insert_right: false,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmfile::parse_cm_file;
    use crate::config::configure_scores;
    use crate::qdb::{calc_qdb_bands, QdbBands};

    /// The QDB-banded alignment must reproduce the dense alignment exactly (the bands contain the
    /// optimal parse), at a fraction of the memory — same score, same per-residue consensus mapping.
    #[test]
    fn banded_alignment_matches_dense() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
        let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
        configure_scores(&mut cm);
        let emap = EmitMap::build(&cm);
        let bands = calc_qdb_bands(&cm, QdbBands::DEFAULT_BETA);

        let cons = ornament_alphabet::read_fasta(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/trna_cons.fa"
        ))
        .expect("read cons")[0]
            .seq
            .clone();
        let dsq = cm.abc.digitize(&cons).expect("digitize");
        let l = dsq.len() - 2;

        let dense = align_glocal(&cm, &dsq, 1, l, &emap);
        let banded = align_glocal_banded(&cm, &dsq, 1, l, &emap, &bands);

        assert_eq!(
            banded.score.to_bits(),
            dense.score.to_bits(),
            "banded score {} != dense {}",
            banded.score,
            dense.score
        );
        assert_eq!(
            banded.residues, dense.residues,
            "banded alignment residues differ from dense"
        );
    }
}
