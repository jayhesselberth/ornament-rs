//! Differential test for glocal CYK alignment + per-residue consensus mapping.
//!
//! Oracle (`cmsearch -g --max --cyk -A`) aligns the RF00005 consensus to all 71 match
//! columns 1:1 (no inserts): the RF line is 71 columns and the aligned sequence has no
//! gaps. The native aligner must reproduce that mapping and the 87.3-bit score.

use easel_rs::{read_fasta, Alphabet};
use infernal_rs::{align_glocal, configure_scores, parse_cm_file, EmitMap};

fn cm() -> infernal_rs::Cm {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/tRNA.cm");
    let mut cm = parse_cm_file(path).expect("parse tRNA.cm");
    configure_scores(&mut cm);
    cm
}

#[test]
fn emit_map_spans_consensus() {
    let cm = cm();
    let emap = EmitMap::build(&cm);
    assert_eq!(emap.clen, 71);
    // For the ROOT (a nonemitter) lpos/rpos are non-inclusive bounds: 0 and clen+1.
    assert_eq!(emap.lpos[0], 0);
    assert_eq!(emap.rpos[0], 72);
    // Every MATP/MATL/MATR match column lies in 1..=clen.
    for nd in 0..cm.nodes {
        if matches!(
            cm.ndtype[nd],
            infernal_rs::model::nd::MATP | infernal_rs::model::nd::MATL
        ) {
            assert!((1..=71).contains(&emap.lpos[nd]), "lpos in range");
        }
        if matches!(
            cm.ndtype[nd],
            infernal_rs::model::nd::MATP | infernal_rs::model::nd::MATR
        ) {
            assert!((1..=71).contains(&emap.rpos[nd]), "rpos in range");
        }
    }
}

#[test]
fn align_consensus_maps_all_columns() {
    let cm = cm();
    let emap = EmitMap::build(&cm);
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_cons.fa");
    let recs = read_fasta(fa).expect("read consensus");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");
    let l = recs[0].seq.len();

    let aln = align_glocal(&cm, &dsq, 1, l, &emap);
    eprintln!("align score = {:.3}", aln.score);
    assert!(
        (aln.score - 87.3).abs() < 0.2,
        "align score {:.3} != oracle 87.3",
        aln.score
    );

    // Oracle aligns all 71 residues as matches to columns 1..71, no inserts.
    assert_eq!(aln.residues.len(), 71);
    let cols: Vec<usize> = aln
        .residues
        .iter()
        .map(|r| r.consensus.expect("every residue is a match column"))
        .collect();
    assert_eq!(cols, (1..=71).collect::<Vec<_>>(), "1:1 column mapping");

    // Residues are in sequence order with strictly increasing consensus columns
    // (tRNA has no pseudoknots), and seq positions are 1..71.
    for (k, r) in aln.residues.iter().enumerate() {
        assert_eq!(r.seq_pos, k + 1);
    }

    // Acceptor-stem / T-loop pairing implies several paired (MP) residues exist.
    assert!(aln.residues.iter().filter(|r| r.paired).count() >= 20);
}

#[test]
fn align_embedded_maps_to_window() {
    // Aligning the embedded hit window (seq 61..131) must yield the same 1:1 column map,
    // now offset into the longer sequence.
    let cm = cm();
    let emap = EmitMap::build(&cm);
    let abc = Alphabet::rna();
    let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trna_embedded.fa");
    let recs = read_fasta(fa).expect("read embedded");
    let dsq = abc.digitize(&recs[0].seq).expect("digitize");

    let aln = align_glocal(&cm, &dsq, 61, 131, &emap);
    assert!((aln.score - 87.3).abs() < 0.2);
    assert_eq!(aln.residues.len(), 71);
    // First residue is at seq position 61 mapped to consensus column 1.
    assert_eq!(aln.residues[0].seq_pos, 61);
    assert_eq!(aln.residues[0].consensus, Some(1));
    assert_eq!(aln.residues[70].seq_pos, 131);
    assert_eq!(aln.residues[70].consensus, Some(71));
}
