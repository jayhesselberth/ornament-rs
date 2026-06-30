//! End-to-end test for `scan --format stockholm` (native alignment output).
//!
//! Drives `infernal::scan_native_aligned` over the multi-tRNA fixture and checks the
//! rendered Stockholm block is well-formed. A second, skip-gated section reruns the same
//! fixture through `cmsearch -A` and asserts parity: identical `#=GC RF`, identical aligned
//! hit rows, and an identical `#=GC SS_cons` (we port Infernal's full-WUSS consensus structure,
//! so the line should match byte-for-byte; the test also compares the decoded base pairs).

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use ornament_core::infernal::{scan_native_aligned, write_stockholm};

const E_VALUE: f64 = 1e-5;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// A parsed Stockholm block: aligned sequence rows (label -> residues), plus the RF and
/// SS_cons GC lines.
#[derive(Default)]
struct Block {
    rows: HashMap<String, String>,
    rf: String,
    ss_cons: String,
}

/// Minimal Stockholm parser for a single block: enough for these fixtures (one block,
/// non-interleaved). Splits each line into `label` + `text`, routing `#=GC RF`/`#=GC SS_cons`
/// to their fields and every other non-comment line to a sequence row.
fn parse_block(text: &str) -> Block {
    let mut b = Block::default();
    for line in text.lines() {
        if line.is_empty() || line == "//" || line == "# STOCKHOLM 1.0" {
            continue;
        }
        if let Some(rest) = line.strip_prefix("#=GC SS_cons") {
            b.ss_cons = rest.trim().to_string();
        } else if let Some(rest) = line.strip_prefix("#=GC RF") {
            b.rf = rest.trim().to_string();
        } else if line.starts_with('#') {
            continue;
        } else {
            let mut it = line.split_whitespace();
            let (Some(label), Some(seq)) = (it.next(), it.next()) else {
                continue;
            };
            b.rows.insert(label.to_string(), seq.to_string());
        }
    }
    b
}

/// Consensus base pairs from a WUSS string as sorted `(i, j)` column indices. Treats every
/// bracket family (`<([{`) as an open and its mirror as a close, nesting via a stack, so two
/// strings that encode the same pairing compare equal regardless of which brackets they use.
fn wuss_pairs(ss: &str) -> Vec<(usize, usize)> {
    let mut stack = Vec::new();
    let mut pairs = Vec::new();
    for (i, c) in ss.chars().enumerate() {
        match c {
            '<' | '(' | '[' | '{' => stack.push(i),
            '>' | ')' | ']' | '}' => {
                let j = stack.pop().expect("balanced WUSS");
                pairs.push((j, i));
            }
            _ => {}
        }
    }
    assert!(stack.is_empty(), "balanced WUSS: {ss}");
    pairs.sort_unstable();
    pairs
}

#[test]
fn stockholm_output_is_well_formed() {
    let msas =
        scan_native_aligned(data("tRNA.cm"), data("trna_multi.fa"), E_VALUE).expect("aligned scan");

    // Single model, three hits (two plus-strand copies + one minus-strand copy in the fixture).
    assert_eq!(msas.len(), 1, "one model block");
    let m = &msas[0];
    assert_eq!(m.rows.len(), 3, "three aligned hits");
    assert_eq!(m.rf.len(), m.clen);
    assert_eq!(m.ss_cons.len(), m.clen);

    let out = write_stockholm(&msas);
    assert!(out.starts_with("# STOCKHOLM 1.0\n"));
    assert!(out.trim_end().ends_with("//"));

    let block = parse_block(&out);
    // Every row, plus RF and SS_cons, render to the same number of columns.
    let width = block.rf.len();
    assert_eq!(block.ss_cons.len(), width);
    for (label, seq) in &block.rows {
        assert_eq!(seq.len(), width, "row {label} width matches RF/SS_cons");
    }
    // SS_cons is balanced and describes a real (non-empty) structure.
    assert!(!wuss_pairs(&block.ss_cons).is_empty());
}

#[test]
fn stockholm_matches_cmsearch_oracle() {
    if Command::new("cmsearch").arg("-h").output().is_err() {
        eprintln!("cmsearch not on PATH; skipping oracle comparison");
        return;
    }

    let cm = data("tRNA.cm");
    let fa = data("trna_multi.fa");
    let aln = std::env::temp_dir().join(format!("ornament_sto_{}.sto", std::process::id()));

    let status = Command::new("cmsearch")
        .args(["-E", "0.00001", "-A"])
        .arg(&aln)
        .arg(&cm)
        .arg(&fa)
        .output()
        .expect("run cmsearch -A");
    assert!(status.status.success(), "cmsearch -A failed");

    let oracle = parse_block(&std::fs::read_to_string(&aln).expect("read oracle .sto"));
    let _ = std::fs::remove_file(&aln);

    let ours = parse_block(&write_stockholm(
        &scan_native_aligned(&cm, &fa, E_VALUE).expect("aligned scan"),
    ));

    // RF (consensus residues) is byte-identical.
    assert_eq!(ours.rf, oracle.rf, "RF line differs from cmsearch");

    // Aligned hit rows match for every shared label.
    for (label, oseq) in &oracle.rows {
        if let Some(ourseq) = ours.rows.get(label) {
            assert_eq!(ourseq, oseq, "aligned residues differ for {label}");
        }
    }
    assert_eq!(
        ours.rows.len(),
        oracle.rows.len(),
        "same number of aligned hits"
    );

    // Full-WUSS SS_cons matches cmsearch byte-for-byte (and, redundantly, the decoded pairs).
    assert_eq!(
        ours.ss_cons, oracle.ss_cons,
        "SS_cons differs from cmsearch"
    );
    assert_eq!(
        wuss_pairs(&ours.ss_cons),
        wuss_pairs(&oracle.ss_cons),
        "consensus base pairs differ from cmsearch"
    );
}
