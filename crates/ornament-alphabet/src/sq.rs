//! Sequences and FASTA input.
//!
//! The digital [`Sequence`] container and its alphabet-aware operations
//! (`digitize`, `reverse_complement`) are ours; FASTA parsing is delegated to
//! [`noodles_fasta`].

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use crate::alphabet::{Alphabet, Dsq};
use crate::AlphabetError;

/// A biosequence with both text and (optionally) digital representations.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: String,
    pub desc: String,
    /// Raw residues as text (no whitespace).
    pub seq: String,
    /// 1-based digital sequence with sentinels at 0 and L+1, if digitized.
    pub dsq: Option<Vec<Dsq>>,
}

impl Sequence {
    /// Residue length L.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Digitize this sequence against `abc`, populating `dsq`.
    pub fn digitize(&mut self, abc: &Alphabet) -> Result<(), AlphabetError> {
        self.dsq = Some(abc.digitize(&self.seq)?);
        Ok(())
    }

    /// Return a reverse-complemented copy (digitizing if needed).
    pub fn reverse_complement(&self, abc: &Alphabet) -> Result<Sequence, AlphabetError> {
        let mut dsq = match &self.dsq {
            Some(d) => d.clone(),
            None => abc.digitize(&self.seq)?,
        };
        abc.revcomp(&mut dsq)?;
        let seq = abc.textize(&dsq);
        Ok(Sequence {
            name: self.name.clone(),
            desc: self.desc.clone(),
            seq,
            dsq: Some(dsq),
        })
    }
}

/// Read all sequences from a FASTA file.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Sequence>, AlphabetError> {
    let file = File::open(path.as_ref())
        .map_err(|e| AlphabetError::Io(format!("{}: {e}", path.as_ref().display())))?;
    read_fasta_from(file)
}

/// Read all FASTA records from any reader, via [`noodles_fasta`].
pub fn read_fasta_from<R: Read>(reader: R) -> Result<Vec<Sequence>, AlphabetError> {
    let mut fasta_reader = noodles_fasta::io::Reader::new(BufReader::new(reader));
    let mut records = Vec::new();

    for result in fasta_reader.records() {
        let rec = result.map_err(|e| AlphabetError::Parse(e.to_string()))?;
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        let desc = rec.description().map(|d| d.to_string()).unwrap_or_default();
        let seq = String::from_utf8_lossy(rec.sequence().as_ref()).into_owned();
        records.push(Sequence {
            name,
            desc,
            seq,
            dsq: None,
        });
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multi_record_fasta() {
        let data = ">seq1 first sequence\nACGU\nACGU\n>seq2\nGGGG\n";
        let recs = read_fasta_from(data.as_bytes()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "seq1");
        assert_eq!(recs[0].desc, "first sequence");
        assert_eq!(recs[0].seq, "ACGUACGU");
        assert_eq!(recs[1].name, "seq2");
        assert_eq!(recs[1].seq, "GGGG");
    }

    #[test]
    fn digitize_and_revcomp_record() {
        let abc = Alphabet::rna();
        let mut recs = read_fasta_from(">s\nAACGU\n".as_bytes()).unwrap();
        recs[0].digitize(&abc).unwrap();
        let rc = recs[0].reverse_complement(&abc).unwrap();
        assert_eq!(rc.seq, "ACGUU");
    }
}
