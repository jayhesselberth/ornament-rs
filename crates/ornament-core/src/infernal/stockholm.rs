//! Stockholm multiple-alignment output for the native scanner (`scan --format stockholm`).
//!
//! Renders the per-hit glocal CYK parses of a model as one Stockholm alignment block —
//! the same shape as `cmsearch -A`: hit rows plus full-WUSS `#=GC SS_cons` and `#=GC RF`. Insert
//! columns are reconciled across all of a model's hits exactly as Infernal's
//! `Parsetrees2Alignment` does: per consensus gap, IL inserts are left-justified and IR
//! inserts right-flushed, and the gap is wide enough for the worst-case hit
//! (`max(IL) + max(IR)`). Each model that produced hits is one `# STOCKHOLM 1.0 … //`
//! block; concatenated blocks (Rfam.seed style) make the multi-model file.

/// One residue in an aligned hit row, in 5'->3' model order.
#[derive(Debug, Clone, Copy)]
pub struct ResCell {
    /// Consensus column (1..=clen) for a match emission; `None` for an insert.
    pub consensus: Option<usize>,
    /// For an insert (`consensus == None`), the consensus column it follows
    /// (IL: `lpos`; IR: `rpos-1`). Ignored for matches.
    pub insert_after: usize,
    /// True for an IR insert (right-flushed within its gap, after IL inserts).
    pub insert_right: bool,
    /// Residue character (uppercase for a match, lowercase for an insert).
    pub ch: char,
}

/// One hit laid out as a row of an alignment.
#[derive(Debug, Clone)]
pub struct AlignedRow {
    /// Row label, `name/start-end`.
    pub label: String,
    /// Residues in 5'->3' model order.
    pub residues: Vec<ResCell>,
}

/// All hits for one model, ready to render as a Stockholm block.
#[derive(Debug, Clone)]
pub struct ModelMsa {
    pub model_name: String,
    pub clen: usize,
    /// Reference (`#=GC RF`) characters, one per consensus column (length `clen`).
    pub rf: Vec<char>,
    /// Consensus structure (`#=GC SS_cons`) characters, one per consensus column.
    pub ss_cons: Vec<char>,
    pub rows: Vec<AlignedRow>,
}

const SS_LABEL: &str = "#=GC SS_cons";
const RF_LABEL: &str = "#=GC RF";

/// Render one or more model alignments as a Stockholm file (one `//`-terminated block per
/// model). Returns an empty string when there are no alignments to write.
pub fn write_stockholm(models: &[ModelMsa]) -> String {
    let mut out = String::new();
    for m in models {
        write_block(&mut out, m);
    }
    out
}

fn write_block(out: &mut String, m: &ModelMsa) {
    let clen = m.clen;

    // Per consensus gap `cpos` in 0..=clen, the worst-case insert widths across all rows.
    // `cpos` is the column an insert follows; cpos 0 holds leading inserts, cpos clen trailing.
    let mut maxil = vec![0usize; clen + 1];
    let mut maxir = vec![0usize; clen + 1];
    for row in &m.rows {
        let mut il = vec![0usize; clen + 1];
        let mut ir = vec![0usize; clen + 1];
        for r in &row.residues {
            if r.consensus.is_some() {
                continue;
            }
            let g = r.insert_after.min(clen);
            if r.insert_right {
                ir[g] += 1;
            } else {
                il[g] += 1;
            }
        }
        for g in 0..=clen {
            maxil[g] = maxil[g].max(il[g]);
            maxir[g] = maxir[g].max(ir[g]);
        }
    }

    out.push_str("# STOCKHOLM 1.0\n");

    let label_w = m
        .rows
        .iter()
        .map(|r| r.label.len())
        .chain([SS_LABEL.len(), RF_LABEL.len()])
        .max()
        .unwrap_or(0);

    // Sequence rows.
    for row in &m.rows {
        let seq = render_row(row, clen, &maxil, &maxir);
        out.push_str(&format!("{:<label_w$}  {}\n", row.label, seq));
    }

    // Annotation rows: SS_cons then RF (insert columns are '.').
    let ss = render_gc(&m.ss_cons, clen, &maxil, &maxir);
    let rf = render_gc(&m.rf, clen, &maxil, &maxir);
    out.push_str(&format!("{:<label_w$}  {}\n", SS_LABEL, ss));
    out.push_str(&format!("{:<label_w$}  {}\n", RF_LABEL, rf));
    out.push_str("//\n");
}

/// Lay a hit row into the reconciled column grid: match residue (uppercase) or `-` per
/// consensus column, IL inserts left-justified and IR inserts right-flushed per gap, `.`
/// padding insert positions this hit doesn't fill.
fn render_row(row: &AlignedRow, clen: usize, maxil: &[usize], maxir: &[usize]) -> String {
    // Bucket this row's residues by consensus column / insert gap.
    let mut match_at = vec![None; clen + 1]; // 1..=clen
    let mut il: Vec<Vec<char>> = vec![Vec::new(); clen + 1];
    let mut ir: Vec<Vec<char>> = vec![Vec::new(); clen + 1];
    for r in &row.residues {
        match r.consensus {
            Some(col) => match_at[col] = Some(r.ch),
            None => {
                let g = r.insert_after.min(clen);
                if r.insert_right {
                    ir[g].push(r.ch);
                } else {
                    il[g].push(r.ch);
                }
            }
        }
    }

    let mut s = String::new();
    for cpos in 0..=clen {
        if cpos >= 1 {
            s.push(match_at[cpos].unwrap_or('-'));
        }
        // IL region: left-justified, '.'-padded.
        for k in 0..maxil[cpos] {
            s.push(il[cpos].get(k).copied().unwrap_or('.'));
        }
        // IR region: right-flushed, '.'-padded on the left.
        let pad = maxir[cpos] - ir[cpos].len();
        for _ in 0..pad {
            s.push('.');
        }
        for &c in &ir[cpos] {
            s.push(c);
        }
    }
    s
}

/// Lay a per-consensus-column annotation (RF / SS_cons) into the same grid: the annotation
/// character at each match column, `.` at every insert position.
fn render_gc(per_col: &[char], clen: usize, maxil: &[usize], maxir: &[usize]) -> String {
    let mut s = String::new();
    for cpos in 0..=clen {
        if cpos >= 1 {
            s.push(per_col.get(cpos - 1).copied().unwrap_or('.'));
        }
        for _ in 0..(maxil[cpos] + maxir[cpos]) {
            s.push('.');
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn match_cell(col: usize, ch: char) -> ResCell {
        ResCell {
            consensus: Some(col),
            insert_after: 0,
            insert_right: false,
            ch,
        }
    }
    fn insert_cell(after: usize, right: bool, ch: char) -> ResCell {
        ResCell {
            consensus: None,
            insert_after: after,
            insert_right: right,
            ch,
        }
    }

    fn lines(block: &str) -> Vec<&str> {
        block.lines().collect()
    }

    #[test]
    fn no_inserts_aligns_to_consensus_columns() {
        let m = ModelMsa {
            model_name: "T".into(),
            clen: 3,
            // RF is emitted verbatim by the writer (uppercasing is the caller's job).
            rf: vec!['A', 'C', 'G'],
            ss_cons: vec!['<', ':', '>'],
            rows: vec![
                AlignedRow {
                    label: "s1/1-3".into(),
                    residues: vec![match_cell(1, 'A'), match_cell(2, 'C'), match_cell(3, 'G')],
                },
                // s2 deletes the middle column.
                AlignedRow {
                    label: "s2/1-2".into(),
                    residues: vec![match_cell(1, 'A'), match_cell(3, 'G')],
                },
            ],
        };
        let out = write_stockholm(&[m]);
        let l = lines(&out);
        assert_eq!(l[0], "# STOCKHOLM 1.0");
        // Three consensus columns, no inserts; deletion is '-'.
        assert!(l[1].ends_with("ACG"), "{}", l[1]);
        assert!(l[2].ends_with("A-G"), "{}", l[2]);
        assert!(l[3].ends_with("<:>"), "{}", l[3]); // SS_cons
        assert!(l[4].ends_with("ACG"), "{}", l[4]); // RF (uppercased by caller-provided chars)
        assert_eq!(*l.last().unwrap(), "//");
    }

    #[test]
    fn insert_widths_reconciled_across_rows() {
        // Both hits insert after column 1; row A has 2 IL, row B has 1 IL -> gap width 2.
        let m = ModelMsa {
            model_name: "T".into(),
            clen: 2,
            rf: vec!['a', 'b'],
            ss_cons: vec![':', ':'],
            rows: vec![
                AlignedRow {
                    label: "a".into(),
                    residues: vec![
                        match_cell(1, 'A'),
                        insert_cell(1, false, 'x'),
                        insert_cell(1, false, 'y'),
                        match_cell(2, 'B'),
                    ],
                },
                AlignedRow {
                    label: "b".into(),
                    residues: vec![
                        match_cell(1, 'A'),
                        insert_cell(1, false, 'z'),
                        match_cell(2, 'B'),
                    ],
                },
            ],
        };
        let out = write_stockholm(&[m]);
        let l = lines(&out);
        // col1 + 2 insert cols + col2
        assert!(l[1].ends_with("AxyB"), "{}", l[1]);
        assert!(l[2].ends_with("Az.B"), "{}", l[2]); // IL left-justified, '.'-padded
                                                     // GC lines mark all insert columns with '.'
        assert!(l[3].ends_with(":..:"), "{}", l[3]);
        assert!(l[4].ends_with("a..b"), "{}", l[4]);
    }

    #[test]
    fn ir_inserts_are_right_flushed_after_il() {
        // One gap after column 1 holds an IL ('x') and, in the other row, an IR ('q').
        // maxil=1, maxir=1 -> gap width 2; IL goes first (left), IR last (right).
        let m = ModelMsa {
            model_name: "T".into(),
            clen: 2,
            rf: vec!['a', 'b'],
            ss_cons: vec![':', ':'],
            rows: vec![
                AlignedRow {
                    label: "il".into(),
                    residues: vec![
                        match_cell(1, 'A'),
                        insert_cell(1, false, 'x'),
                        match_cell(2, 'B'),
                    ],
                },
                AlignedRow {
                    label: "ir".into(),
                    residues: vec![
                        match_cell(1, 'A'),
                        insert_cell(1, true, 'q'),
                        match_cell(2, 'B'),
                    ],
                },
            ],
        };
        let out = write_stockholm(&[m]);
        let l = lines(&out);
        // IL row: IL filled, IR region '.' -> "Ax.B"
        assert!(l[1].ends_with("Ax.B"), "{}", l[1]);
        // IR row: IL region '.', IR right-flushed -> "A.qB"
        assert!(l[2].ends_with("A.qB"), "{}", l[2]);
    }

    #[test]
    fn multiple_models_emit_multiple_blocks() {
        let mk = |name: &str| ModelMsa {
            model_name: name.into(),
            clen: 1,
            rf: vec!['a'],
            ss_cons: vec![':'],
            rows: vec![AlignedRow {
                label: format!("{name}/1-1"),
                residues: vec![match_cell(1, 'A')],
            }],
        };
        let out = write_stockholm(&[mk("m1"), mk("m2")]);
        assert_eq!(out.matches("# STOCKHOLM 1.0").count(), 2);
        assert_eq!(out.matches("\n//\n").count(), 2);
    }
}
