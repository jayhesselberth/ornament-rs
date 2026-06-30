//! Native Sprinzl-position assignment via the reference covariance model.
//!
//! This replaces the old broken 1:1 fallback in `analysis::compatibility`. A tRNA sequence
//! is aligned to a Sprinzl-numbered reference CM (built from clover ground truth; see
//! `scripts/build-sprinzl-cm.py`) with the native `ornament-scfg` glocal aligner, and each
//! matched residue is mapped to its consensus column → Sprinzl label.
//!
//! The reference model and its column→label table are embedded so annotation needs no
//! external files or the `cmsearch`/`cmalign` binaries.

use std::collections::HashMap;
use std::sync::OnceLock;

use ornament_alphabet::Alphabet;
use ornament_scfg::{align_glocal, configure_scores, parse_cm_str, Cm, EmitMap};

use super::types::SprinzlPosition;

const SPRINZL_CM: &str = include_str!("data/sprinzl_euk.cm");
const GCOL_LABEL: &str = include_str!("data/sprinzl_gcol_label.tsv");

struct SprinzlModel {
    cm: Cm,
    emap: EmitMap,
    /// Consensus column (1..clen) → Sprinzl label.
    col_label: HashMap<usize, String>,
    abc: Alphabet,
}

fn model() -> &'static SprinzlModel {
    static M: OnceLock<SprinzlModel> = OnceLock::new();
    M.get_or_init(|| {
        let mut cm = parse_cm_str(SPRINZL_CM).expect("embedded Sprinzl CM parses");
        configure_scores(&mut cm);
        let emap = EmitMap::build(&cm);
        let mut col_label = HashMap::new();
        for line in GCOL_LABEL.lines().skip(1) {
            let mut it = line.split('\t');
            if let (Some(c), Some(l)) = (it.next(), it.next()) {
                if let Ok(c) = c.parse::<usize>() {
                    col_label.insert(c, l.to_string());
                }
            }
        }
        SprinzlModel {
            cm,
            emap,
            col_label,
            abc: Alphabet::rna(),
        }
    })
}

/// Align a tRNA sequence to the Sprinzl reference model, returning a map from Sprinzl
/// position to the 0-based index into `sequence`. Insert residues (no consensus column)
/// are omitted. Returns an empty map if the sequence can't be digitized.
pub fn align_to_sprinzl(sequence: &str) -> HashMap<SprinzlPosition, usize> {
    let mut out = HashMap::new();
    if sequence.is_empty() {
        return out;
    }
    let m = model();
    let dsq = match m.abc.digitize(sequence) {
        Ok(d) => d,
        Err(_) => return out, // non-standard residues: caller falls back
    };
    let aln = align_glocal(&m.cm, &dsq, 1, sequence.len(), &m.emap);
    for r in &aln.residues {
        if let Some(col) = r.consensus {
            if let Some(lab) = m.col_label.get(&col) {
                out.insert(SprinzlPosition(lab.clone()), r.seq_pos - 1);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligns_yeast_phe_to_anticodon() {
        // S. cerevisiae tRNA-Phe (standard length). Anticodon (34-36) = GAA.
        let seq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";
        let m = align_to_sprinzl(seq);
        // Core landmarks present and mapping into the sequence.
        let p34 = m.get(&SprinzlPosition("34".into())).copied();
        let p35 = m.get(&SprinzlPosition("35".into())).copied();
        let p36 = m.get(&SprinzlPosition("36".into())).copied();
        assert!(
            p34.is_some() && p35.is_some() && p36.is_some(),
            "anticodon mapped"
        );
        // Anticodon residues are consecutive.
        assert_eq!(p35, p34.map(|x| x + 1));
        assert_eq!(p36, p35.map(|x| x + 1));
        // Position 1 maps to the first residue.
        assert_eq!(m.get(&SprinzlPosition("1".into())).copied(), Some(0));
    }
}
