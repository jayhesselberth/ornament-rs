//! Generalization test: the native engine on the tRNAscan-SE eukaryotic model
//! (`TRNAinf-euk.cm`, INFERNAL1/a 1.1.1, CLEN=90, W=628) — a different, larger model than
//! RF00005. Oracle (`cmsearch -g --max [--cyk]`): the 90-nt consensus scores 114.0 (CYK) /
//! 114.2 (Inside) over seq 1..90.

use easel_rs::{read_fasta, Alphabet};
use infernal_rs::{
    align_glocal, configure_scores, cyk_scan_glocal, inside_scan_glocal, parse_cm_file, EmitMap,
};

fn setup() -> (infernal_rs::Cm, Vec<easel_rs::Dsq>, usize) {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/TRNAinf-euk.cm");
    let mut cm = parse_cm_file(path).expect("parse TRNAinf-euk.cm");
    assert_eq!(cm.clen, 90, "euk model consensus length");
    assert_eq!(cm.m, 289, "euk model states");
    configure_scores(&mut cm);
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/euk_cons.fa");
    let recs = read_fasta(fa).expect("read euk consensus");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");
    let l = recs[0].seq.len();
    (cm, dsq, l)
}

#[test]
fn euk_model_cyk_and_inside_match_oracle() {
    let (cm, dsq, l) = setup();
    assert_eq!(l, 90);

    let cyk = cyk_scan_glocal(&cm, &dsq, cm.w as usize).best.expect("cyk hit");
    eprintln!("euk CYK: score={:.3} {}..{}", cyk.score, cyk.i, cyk.j);
    assert_eq!((cyk.i, cyk.j), (1, 90));
    assert!((cyk.score - 114.0).abs() < 0.2, "CYK {:.3} != 114.0", cyk.score);

    let ins = inside_scan_glocal(&cm, &dsq, cm.w as usize).best.expect("inside hit");
    eprintln!("euk Inside: score={:.3}", ins.score);
    assert!((ins.score - 114.2).abs() < 0.25, "Inside {:.3} != 114.2", ins.score);
    assert!(ins.score >= cyk.score - 1e-3);
}

#[test]
fn euk_model_alignment_matches_oracle_structure() {
    // The oracle (cmsearch -g --max --cyk -A) aligns this consensus with exactly one
    // insert (a lowercase residue) and one deleted match column (a '-' gap) — it is NOT a
    // clean 1:1 map. The native aligner must reproduce that exact structure.
    let (cm, dsq, l) = setup();
    let emap = EmitMap::build(&cm);
    assert_eq!(emap.clen, 90);
    let aln = align_glocal(&cm, &dsq, 1, l, &emap);
    assert!((aln.score - 114.0).abs() < 0.2, "align score {:.3}", aln.score);

    assert_eq!(aln.residues.len(), 90, "all 90 consensus residues placed");

    // Exactly one residue is an insert (no consensus column).
    let inserts = aln.residues.iter().filter(|r| r.consensus.is_none()).count();
    assert_eq!(inserts, 1, "oracle has exactly one inserted residue");

    // Match columns are strictly increasing, within 1..=90, with exactly one deleted.
    let cols: Vec<usize> = aln.residues.iter().filter_map(|r| r.consensus).collect();
    assert_eq!(cols.len(), 89, "89 occupied match columns (one deleted)");
    assert!(cols.windows(2).all(|w| w[0] < w[1]), "columns strictly increasing");
    assert!(cols.iter().all(|&c| (1..=90).contains(&c)), "columns in range");
    let missing: Vec<usize> = (1..=90).filter(|c| !cols.contains(c)).collect();
    assert_eq!(missing.len(), 1, "exactly one deleted consensus column");
}
