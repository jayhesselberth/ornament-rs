//! Sequences and FASTA input (the `esl_sq` / `esl_sqio` subset we need).

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::alphabet::{Alphabet, Dsq};
use crate::EaselError;

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
    pub fn digitize(&mut self, abc: &Alphabet) -> Result<(), EaselError> {
        self.dsq = Some(abc.digitize(&self.seq)?);
        Ok(())
    }

    /// Return a reverse-complemented copy (digitizing if needed).
    pub fn reverse_complement(&self, abc: &Alphabet) -> Result<Sequence, EaselError> {
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
pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Sequence>, EaselError> {
    let file = File::open(path.as_ref())
        .map_err(|e| EaselError::Io(format!("{}: {e}", path.as_ref().display())))?;
    read_fasta_from(BufReader::new(file))
}

/// Read all FASTA records from any reader.
pub fn read_fasta_from<R: Read>(reader: R) -> Result<Vec<Sequence>, EaselError> {
    let buf = BufReader::new(reader);
    let mut records = Vec::new();
    let mut cur: Option<Sequence> = None;

    for line in buf.lines() {
        let line = line.map_err(|e| EaselError::Io(e.to_string()))?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some(prev) = cur.take() {
                records.push(prev);
            }
            let header = header.trim_start();
            let (name, desc) = match header.split_once(char::is_whitespace) {
                Some((n, d)) => (n.to_string(), d.trim().to_string()),
                None => (header.to_string(), String::new()),
            };
            cur = Some(Sequence {
                name,
                desc,
                seq: String::new(),
                dsq: None,
            });
        } else {
            match cur.as_mut() {
                Some(rec) => rec.seq.push_str(trimmed.trim()),
                None => return Err(EaselError::Io("FASTA data before first '>' header".into())),
            }
        }
    }
    if let Some(prev) = cur.take() {
        records.push(prev);
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
