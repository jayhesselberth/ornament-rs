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

#[test]
fn assigns_sprinzl_labels_matching_clover_truth() {
    let mut cm = parse_cm_file(data("sprinzl_euk.cm")).expect("parse");
    configure_scores(&mut cm);
    let emap = EmitMap::build(&cm);
    let labels = column_labels();
    let abc = Alphabet::rna();

    let table = std::fs::read_to_string(data("sprinzl_test_trnas.tsv")).unwrap();
    let mut checked = 0;
    for line in table.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let (id, seq, truth) = (cols[0], cols[1], cols[2]);
        let truth: Vec<&str> = truth.split(',').collect();
        assert_eq!(seq.len(), truth.len());

        let dsq = abc.digitize(seq).expect("digitize");
        let aln = align_glocal(&cm, &dsq, 1, seq.len(), &emap);

        // Map each residue's 1-based seq position -> assigned Sprinzl label (match cols).
        let mut assigned: HashMap<usize, String> = HashMap::new();
        for r in &aln.residues {
            if let Some(col) = r.consensus {
                if let Some(lab) = labels.get(&col) {
                    assigned.insert(r.seq_pos, lab.clone());
                }
            }
        }

        // Compare against ground truth at each position. These are standard-length tRNAs,
        // so every residue should be a match column with the correct Sprinzl label.
        let mut correct = 0;
        for (i, want) in truth.iter().enumerate() {
            let pos = i + 1;
            if assigned.get(&pos).map(String::as_str) == Some(*want) {
                correct += 1;
            }
        }
        let frac = correct as f64 / truth.len() as f64;
        eprintln!("{id}: {correct}/{} Sprinzl labels correct ({:.1}%)", truth.len(), frac * 100.0);
        assert!(
            frac >= 0.95,
            "{id}: only {correct}/{} Sprinzl labels match clover truth",
            truth.len()
        );
        checked += 1;
    }
    assert!(checked >= 2, "expected at least 2 test tRNAs");
}
