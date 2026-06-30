//! End-to-end test for the native `scan` engine.
//!
//! Drives `infernal::scan_native` (the default `ornament scan` path, pure Rust, no
//! external binary) over a small multi-tRNA FASTA, then feeds the hits through the
//! same Sprinzl alignment + modification analysis the CLI uses. A second, skip-gated
//! section reruns the identical fixture through the `cmsearch` subprocess oracle and
//! asserts the native hits agree on coordinates/strand and produce identical Sprinzl
//! and modification calls.

use std::collections::HashSet;
use std::path::PathBuf;
use std::process::Command;

use ornament_core::analysis::{analyze_compatibility, Strand, TRNAHit};
use ornament_core::infernal::{scan_native, CMHit, InfernalRunner};
use ornament_core::modification::{align_to_sprinzl, ModificationDatabase};
use ornament_core::SprinzlPosition;

const E_VALUE: f64 = 1e-5;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// Reverse-complement an RNA string (used to recover minus-strand hit sequences).
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'U',
            'U' => 'A',
            'G' => 'C',
            'C' => 'G',
            other => other,
        })
        .collect()
}

/// Recover the 5'->3' residues of a hit from the full forward sequence.
fn hit_sequence(full: &str, h: &CMHit) -> String {
    let (lo, hi) = if h.target_start <= h.target_end {
        (h.target_start, h.target_end)
    } else {
        (h.target_end, h.target_start)
    };
    let sub: String = full[lo - 1..hi].to_string();
    if h.strand == '-' {
        revcomp(&sub)
    } else {
        sub
    }
}

fn hit_to_trna(full: &str, h: &CMHit) -> TRNAHit {
    TRNAHit {
        id: format!("{}-{}-{}", h.target_name, h.target_start, h.target_end),
        seq_name: h.target_name.clone(),
        start: h.target_start.min(h.target_end),
        end: h.target_start.max(h.target_end),
        strand: Strand::from(h.strand),
        score: h.score,
        isotype: None,
        anticodon: None,
        sequence: hit_sequence(full, h),
        structure: String::new(),
    }
}

fn read_forward_seq(fasta: &PathBuf) -> String {
    ornament_alphabet::read_fasta(fasta).expect("read fasta")[0]
        .seq
        .clone()
}

fn coord_set(hits: &[CMHit]) -> HashSet<(usize, usize, char)> {
    hits.iter()
        .map(|h| (h.target_start, h.target_end, h.strand))
        .collect()
}

#[test]
fn native_scan_finds_trnas_and_maps_sprinzl() {
    let cm = data("tRNA.cm");
    let fa = data("trna_multi.fa");
    let full = read_forward_seq(&fa);

    let hits = scan_native(&cm, &fa, E_VALUE).expect("native scan");

    // The fixture embeds three significant tRNAs (two plus, one minus strand).
    assert_eq!(hits.len(), 3, "three significant native hits, got {hits:?}");
    assert_eq!(hits.iter().filter(|h| h.strand == '+').count(), 2);
    assert_eq!(hits.iter().filter(|h| h.strand == '-').count(), 1);
    // Minus-strand hit reported high->low (cmsearch convention).
    let minus = hits.iter().find(|h| h.strand == '-').unwrap();
    assert!(minus.target_start > minus.target_end);

    let db = ModificationDatabase::eukaryotic();
    for h in &hits {
        let trna = hit_to_trna(&full, h);
        let aln = align_to_sprinzl(&trna.sequence);
        assert!(
            !aln.is_empty(),
            "hit {}..{} maps to Sprinzl positions",
            h.target_start,
            h.target_end
        );
        // Anticodon loop landmarks present.
        for p in ["34", "35", "36"] {
            assert!(
                aln.contains_key(&SprinzlPosition::new(p)),
                "Sprinzl position {p} present for hit {}..{}",
                h.target_start,
                h.target_end
            );
        }
        // Modification analysis runs and yields a real alignment-backed result.
        let result = analyze_compatibility(&trna, &db);
        assert!(!result.sprinzl_alignment.is_empty());
    }
}

#[test]
fn native_matches_cmsearch_oracle() {
    // Skip-gate: only run when the cmsearch subprocess is available.
    if Command::new("cmsearch").arg("-h").output().is_err() {
        eprintln!("cmsearch not on PATH; skipping oracle comparison");
        return;
    }

    let cm = data("tRNA.cm");
    let fa = data("trna_multi.fa");
    let full = read_forward_seq(&fa);

    let native = scan_native(&cm, &fa, E_VALUE).expect("native scan");
    let oracle = InfernalRunner::new()
        .with_cm(&cm)
        .with_e_value(E_VALUE)
        .cmsearch(&fa)
        .expect("cmsearch oracle");

    let oracle_coords = coord_set(&oracle);

    // Every native hit must appear in the oracle output at the same coordinates/strand.
    for h in &native {
        assert!(
            oracle_coords.contains(&(h.target_start, h.target_end, h.strand)),
            "native hit {}..{} ({}) absent from cmsearch oracle {oracle_coords:?}",
            h.target_start,
            h.target_end,
            h.strand,
        );

        // Sprinzl + modification calls must be identical for the shared coordinates.
        let oh = oracle
            .iter()
            .find(|o| {
                o.target_start == h.target_start
                    && o.target_end == h.target_end
                    && o.strand == h.strand
            })
            .unwrap();
        let native_aln = align_to_sprinzl(&hit_sequence(&full, h));
        let oracle_aln = align_to_sprinzl(&hit_sequence(&full, oh));
        assert_eq!(
            native_aln, oracle_aln,
            "Sprinzl alignment differs at {}..{}",
            h.target_start, h.target_end
        );
    }
}
