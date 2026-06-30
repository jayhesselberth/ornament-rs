//! Consensus emit map (`CMEmitMap_t`): maps model nodes to consensus columns.
//!
//! Faithful port of `CreateEmitMap` (`src/display.c`). For a MATL/MATP node the left match
//! column is `lpos[nd]`; for a MATR/MATP node the right match column is `rpos[nd]`. This is
//! what turns a parse into per-residue consensus (Sprinzl) positions.

use crate::model::{nd, Cm};

/// Per-node consensus bounds.
#[derive(Debug, Clone)]
pub struct EmitMap {
    /// Left consensus bound of the subtree under node `nd`.
    pub lpos: Vec<i32>,
    /// Right consensus bound of the subtree under node `nd`.
    pub rpos: Vec<i32>,
    /// EL inserts come after this consensus position.
    pub epos: Vec<i32>,
    pub clen: usize,
}

impl EmitMap {
    /// Build the emit map by the same left/right depth-first traversal Infernal uses.
    pub fn build(cm: &Cm) -> EmitMap {
        let nodes = cm.nodes;
        let mut lpos = vec![-1i32; nodes];
        let mut rpos = vec![-1i32; nodes];
        let mut epos = vec![-1i32; nodes];

        // Stack of (on_right, nd). on_right=0 -> visiting left side, 1 -> right side.
        let mut pda: Vec<(u8, usize)> = Vec::new();
        let mut cpos: i32 = 0;
        pda.push((0, 0));

        while let Some((on_right, node)) = pda.pop() {
            if on_right == 1 {
                rpos[node] = cpos + 1;
                if cm.ndtype[node] == nd::MATP || cm.ndtype[node] == nd::MATR {
                    cpos += 1;
                }
            } else {
                if cm.ndtype[node] == nd::MATP || cm.ndtype[node] == nd::MATL {
                    cpos += 1;
                }
                lpos[node] = cpos;

                if cm.ndtype[node] == nd::BIF {
                    // The bifurcation's B state holds its BEGL/BEGR children.
                    let b_state = cm.nodemap[node];
                    let left_child_nd = cm.ndidx[cm.cfirst[b_state] as usize];
                    let right_child_nd = cm.ndidx[cm.cnum[b_state] as usize];
                    // Revisit this node on the right after both children.
                    pda.push((1, node));
                    // Right child, then left child (LIFO => left processed first).
                    pda.push((0, right_child_nd));
                    pda.push((0, left_child_nd));
                } else {
                    pda.push((1, node)); // revisit on the right
                    if cm.ndtype[node] != nd::END {
                        pda.push((0, node + 1)); // single child node is nd+1
                    }
                }
            }
        }

        // epos: the consensus position an EL transition from node `nd` follows.
        let mut c: i32 = 0;
        for node in (0..nodes).rev() {
            if cm.ndtype[node] == nd::END {
                c = lpos[node];
            } else if cm.ndtype[node] == nd::BIF {
                let b_state = cm.nodemap[node];
                let right_child_nd = cm.ndidx[cm.cnum[b_state] as usize];
                c = epos[right_child_nd];
            }
            epos[node] = c;
        }

        EmitMap {
            lpos,
            rpos,
            epos,
            clen: cm.clen,
        }
    }
}
