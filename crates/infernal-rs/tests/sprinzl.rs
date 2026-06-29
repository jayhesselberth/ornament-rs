//! End-to-end validation of Sprinzl-position assignment via a purpose-built reference CM.
//!
//! Background: tRNAscan-SE's `*_sprinzl.cm` reference models were never released (confirmed
//! by the authors in UCSC-LoweLab/tRNAscan-SE#14). We rebuild one from the clover project's
//! Sprinzl ground truth (`sacCer_global_coords.tsv.gz`): a 115-column all-match model
//! (`cmbuild --hand`) whose consensus column *i* is exactly Sprinzl label *i*.
//!
//! This test aligns real yeast tRNAs to that reference with the native glocal aligner and
//! checks the assigned Sprinzl labels against clover's ground truth.

use std::collections::HashMap;

use easel_rs::Alphabet;
use infernal_rs::{align_glocal, configure_scores, parse_cm_file, EmitMap};

fn data(name: &str) -> String {
    format!("{}/tests/data/{}", env!("CARGO_MANIFEST_DIR"), name)
}

/// consensus column (1..clen) -> Sprinzl label. With the all-match `--hand` build, consensus
/// column i corresponds to global/seed column i, whose label is in the gcol table.
fn column_labels() -> HashMap<usize, String> {
    let txt = std::fs::read_to_string(data("sprinzl_gcol_label.tsv")).unwrap();
    let mut m = HashMap::new();
    for line in txt.lines().skip(1) {
        let mut it = line.split('\t');
        if let (Some(c), Some(l)) = (it.next(), it.next()) {
            if let Ok(c) = c.parse::<usize>() {
                m.insert(c, l.to_string());
            }
        }
    }
    m
}

#[test]
fn reference_cm_parses_and_is_all_match() {
    let mut cm = parse_cm_file(data("sprinzl_euk.cm")).expect("parse sprinzl_euk.cm");
    assert_eq!(cm.clen, 115, "all-match reference has 115 consensus columns");
    configure_scores(&mut cm);
    // 115 global columns label the full Sprinzl set incl. variable positions; global
    // column index != Sprinzl number (variable columns 17/17a/20a/e* shift it). Check the
    // table covers the canonical landmarks somewhere.
    let labels = column_labels();
    assert_eq!(labels.get(&1).map(String::as_str), Some("1"));
    let set: std::collections::HashSet<&str> = labels.values().map(String::as_str).collect();
    // Anticodon (34-36), a D-loop insertion (20a), variable-loop (e11), discriminator (73)
    // and CCA (76) all appear in the yeast-derived column set.
    for landmark in ["34", "35", "36", "20a", "e11", "73", "76"] {
        assert!(set.contains(landmark), "Sprinzl landmark {landmark} present");
    }
}

/// Assign Sprinzl labels to a sequence via the native aligner: seq position -> label.
fn assign(
    cm: &infernal_rs::Cm,
    emap: &EmitMap,
    labels: &HashMap<usize, String>,
    abc: &Alphabet,
    seq: &str,
) -> HashMap<usize, String> {
    let dsq = abc.digitize(seq).expect("digitize");
    let aln = align_glocal(cm, &dsq, 1, seq.len(), emap);
    let mut out = HashMap::new();
    for r in &aln.residues {
        if let Some(col) = r.consensus {
            if let Some(lab) = labels.get(&col) {
                out.insert(r.seq_pos, lab.clone());
            }
        }
    }
    out
}

fn is_variable(label: &str) -> bool {
    label.starts_with('e') || label.ends_with('a') || label.ends_with('b') || label.is_empty()
}

#[test]
fn standard_trnas_match_clover_exactly() {
    // Standard-length tRNAs (no long variable arm, canonical D-loop) must reproduce clover's
    // Sprinzl numbering exactly — the strong correctness check for the whole pipeline.
    let mut cm = parse_cm_file(data("sprinzl_euk.cm")).expect("parse");
    configure_scores(&mut cm);
    let emap = EmitMap::build(&cm);
    let labels = column_labels();
    let abc = Alphabet::rna();

    let table = std::fs::read_to_string(data("sprinzl_test_trnas.tsv")).unwrap();
    let mut checked = 0;
    for line in table.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        let (id, seq, truth) = (c[0], c[1], c[2].split(',').collect::<Vec<_>>());
        // Only the standard tRNAs here (no variable-region labels in truth).
        if truth.iter().any(|l| is_variable(l)) {
            continue;
        }
        let assigned = assign(&cm, &emap, &labels, &abc, seq);
        let correct = truth
            .iter()
            .enumerate()
            .filter(|(i, w)| assigned.get(&(i + 1)).map(String::as_str) == Some(**w))
            .count();
        eprintln!("{id}: {correct}/{} (standard)", truth.len());
        assert_eq!(correct, truth.len(), "{id}: standard tRNA must match clover exactly");
        checked += 1;
    }
    assert!(checked >= 2);
}

#[test]
fn variable_arm_trnas_correct_at_conserved_positions() {
    // Long-variable-arm tRNAs (Leu, Ser) have variable D-loop / variable-arm regions where a
    // structural CM alignment legitimately differs from clover's per-tRNA convention (e.g.
    // which D-stem position is "deleted"). But the conserved, modification-relevant positions
    // must still be assigned correctly. Validate those by label, not raw position.
    let mut cm = parse_cm_file(data("sprinzl_euk.cm")).expect("parse");
    configure_scores(&mut cm);
    let emap = EmitMap::build(&cm);
    let labels = column_labels();
    let abc = Alphabet::rna();

    // Modification-relevant / conserved landmarks (anticodon loop, T-loop, acceptor stem).
    let conserved = [
        "1", "2", "26", "27", "31", "32", "33", "34", "35", "36", "37", "38", "39", "48", "49",
        "53", "54", "55", "57", "58", "61", "73",
    ];

    let table = std::fs::read_to_string(data("sprinzl_test_trnas.tsv")).unwrap();
    let mut checked = 0;
    for line in table.lines().skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        let (id, seq, truth) = (c[0], c[1], c[2].split(',').collect::<Vec<_>>());
        if !truth.iter().any(|l| is_variable(l)) {
            continue; // only the variable-arm tRNAs
        }
        let assigned = assign(&cm, &emap, &labels, &abc, seq);
        // clover seq position for each label.
        let truth_pos: HashMap<&str, usize> =
            truth.iter().enumerate().map(|(i, l)| (*l, i + 1)).collect();

        let mut ok = 0;
        let mut tot = 0;
        for lab in conserved {
            if let Some(&pos) = truth_pos.get(lab) {
                tot += 1;
                if assigned.get(&pos).map(String::as_str) == Some(lab) {
                    ok += 1;
                }
            }
        }
        eprintln!("{id}: {ok}/{tot} conserved landmarks correct");
        assert_eq!(ok, tot, "{id}: all conserved landmarks must match clover");
        checked += 1;
    }
    assert!(checked >= 1, "expected at least one variable-arm tRNA");
}
