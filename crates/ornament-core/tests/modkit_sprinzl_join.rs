//! End-to-end integration test for the modkit BedMethyl -> Sprinzl join.
//!
//! Builds a synthetic bedMethyl file against a yeast tRNA-Phe hit and checks that modkit
//! calls (a) parse, (b) map to the correct Sprinzl positions through the native reference
//! alignment, and (c) drive the odd-tRNA / compatibility verdict.

use std::collections::HashMap;

use ornament_core::analysis::{
    analyze_compatibility, analyze_compatibility_with_mods, MeasuredVerdict, Strand, TRNAHit,
};
use ornament_core::integration::modkit::{join_to_sprinzl, parse_bedmethyl};
use ornament_core::modification::{align_to_sprinzl, ModificationDatabase};
use ornament_core::SprinzlPosition;

// S. cerevisiae nuc-tRNA-Phe-GAA-1-1 (standard length, 73 nt) from the Sprinzl test set.
const PHE_SEQ: &str = "GCGGACUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAGUUCGCA";

fn phe_hit() -> TRNAHit {
    TRNAHit {
        id: "phe-test".to_string(),
        seq_name: "chrT".to_string(),
        // 1-based inclusive cmsearch-style coords: seq[0] sits at genomic position 1.
        start: 1,
        end: PHE_SEQ.len(),
        strand: Strand::Plus,
        score: 90.0,
        isotype: Some("Phe".to_string()),
        anticodon: Some("GAA".to_string()),
        sequence: PHE_SEQ.to_string(),
        structure: String::new(),
    }
}

/// Build one extended (18-column) modkit bedMethyl line. On the plus strand with the hit
/// anchored at genomic position 1, a residue at sequence index `i` has 0-based BED start `i`.
fn bedmethyl_line(seq_idx: usize, mod_code: &str, coverage: u32, percent: f64) -> String {
    format!(
        "chrT\t{start}\t{end}\t{code}\t{cov}\t+\t{start}\t{end}\t0,0,0\t{cov}\t{pct:.2}\t{cov}\t0\t0\t0\t0\t0\t0",
        start = seq_idx,
        end = seq_idx + 1,
        code = mod_code,
        cov = coverage,
        pct = percent,
    )
}

#[test]
fn modkit_calls_join_to_sprinzl_and_drive_oddness() {
    let hit = phe_hit();
    let db = ModificationDatabase::eukaryotic();
    let aln = align_to_sprinzl(&hit.sequence); // Sprinzl label -> 0-based seq index

    // -- Pick two target positions from the real alignment ------------------------------
    // (a) Position 55: universally pseudouridylated; base must be U here.
    let pos55 = SprinzlPosition::new("55");
    let idx55 = *aln.get(&pos55).expect("Phe aligns position 55");
    assert_eq!(
        hit.sequence.chars().nth(idx55),
        Some('U'),
        "position 55 is uridine (pseudouridine-compatible)"
    );

    // (b) A position with NO MODOMICS expectation whose base is G, so a measured m1A call
    //     (which requires A) is impossible chemistry that sequence-only analysis can't see.
    let iso = ornament_core::modification::Isotype::new("Phe");
    let (g_pos, g_idx) = {
        let mut candidates: Vec<(SprinzlPosition, usize)> = aln
            .iter()
            .filter(|(pos, &idx)| {
                hit.sequence.chars().nth(idx) == Some('G')
                    && db.get_expectations_for_isotype(pos, &iso).is_empty()
            })
            .map(|(pos, &idx)| (pos.clone(), idx))
            .collect();
        candidates.sort_by_key(|(_, idx)| *idx);
        candidates
            .into_iter()
            .next()
            .expect("a G with no expectation")
    };

    // -- Synthetic bedMethyl file --------------------------------------------------------
    let consistent_line = bedmethyl_line(idx55, "17802", 60, 92.0); // pseudouridine (ChEBI)
    let incompatible_line = bedmethyl_line(g_idx, "m1A", 50, 88.0); // m1A on a G: impossible
    let content = format!("{consistent_line}\n{incompatible_line}\n");

    let records = parse_bedmethyl(&content);
    assert_eq!(records.len(), 2, "both bedMethyl rows parse");
    assert_eq!(records[0].mod_code, "17802");
    assert_eq!(records[0].coverage, 60);
    assert!((records[0].mod_frequency - 92.0).abs() < 1e-9);

    // -- Join onto Sprinzl positions -----------------------------------------------------
    let observed = join_to_sprinzl(&hit, &records);
    assert_eq!(observed.len(), 2, "both calls map to a Sprinzl position");
    assert_eq!(
        observed.get(&pos55).map(|c| c.mod_code.as_str()),
        Some("17802"),
        "Psi call lands on position 55"
    );
    assert_eq!(
        observed.get(&g_pos).map(|c| c.mod_code.as_str()),
        Some("m1A"),
        "m1A call lands on the chosen G position"
    );

    // -- Feed into analysis --------------------------------------------------------------
    let base_only = analyze_compatibility(&hit, &db);

    // A single consistent call must not change the oddness verdict.
    let mut just_consistent: HashMap<SprinzlPosition, _> = HashMap::new();
    just_consistent.insert(pos55.clone(), observed.get(&pos55).unwrap().clone());
    let consistent_result = analyze_compatibility_with_mods(&hit, &db, &just_consistent);
    assert_eq!(
        consistent_result.is_odd, base_only.is_odd,
        "a consistent measured mod doesn't change oddness"
    );
    let v55 = consistent_result
        .measured
        .iter()
        .find(|m| m.position == pos55)
        .expect("position 55 graded");
    assert_eq!(v55.verdict, MeasuredVerdict::Consistent);

    // Adding the chemically-impossible m1A call must flag the tRNA as odd.
    let full = analyze_compatibility_with_mods(&hit, &db, &observed);
    assert!(
        full.is_odd,
        "an incompatible measured mod makes the tRNA odd"
    );
    let vg = full
        .measured
        .iter()
        .find(|m| m.position == g_pos)
        .expect("G position graded");
    assert_eq!(
        vg.verdict,
        MeasuredVerdict::Incompatible,
        "m1A on a G is graded incompatible"
    );
    assert_eq!(vg.observed_base, ornament_core::RnaBase::G);
}
