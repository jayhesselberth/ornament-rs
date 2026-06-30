//! Biosequence alphabets with digital-sequence encoding.
//!
//! This is a faithful port of the Easel alphabet machinery (`esl_alphabet`),
//! **extended** with first-class modified-base symbols for modification-aware
//! covariance-model search.
//!
//! ## Layout (faithful to Easel)
//!
//! A nucleic alphabet stores `K` canonical residues (ACGU, `K=4`) and `Kp` total
//! codes. Codes are laid out as:
//!
//! ```text
//!   0 .. K-1     canonical residues          (ACGU)
//!   K            gap                          ('-')
//!   K+1 .. ...   standard IUPAC degeneracies  (RYMKSWHBVD)
//!   ...          registered modified symbols  (Ψ, D, I, m1A, ...)   <-- the extension
//!   Kp-3         "any" / unknown              ('N')
//!   Kp-2         nonresidue                   ('*')
//!   Kp-1         missing                      ('~')
//! ```
//!
//! The `degen[code][0..K]` matrix records which canonical residues a code marginalizes
//! over. A modified symbol's degeneracy vector is the indicator of its **parent**
//! canonical base, so `favg_score(modified) == sc[parent]` exactly — homology scoring is
//! numerically identical to the stock 4-symbol alphabet, while the modified identity is
//! retained for downstream annotation. Phase B (discovery) can later promote modified
//! symbols to their own emission distributions; nothing here precludes that.

use std::collections::HashMap;

use crate::AlphabetError;

/// A digital sequence residue code (`ESL_DSQ`).
pub type Dsq = u8;

/// Sentinel byte at positions 0 and L+1 of a digital sequence.
pub const SENTINEL: Dsq = 255;
/// Code for an input symbol that is unmapped/unexpected.
pub const ILLEGAL: Dsq = 254;
/// `1 - e^x ~ -x` cutoff, matching Easel's `eslSMALLX1`.
pub const SMALLX1: f64 = 5e-9;

/// Standard alphabet kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlphabetType {
    Rna,
    Dna,
    Custom,
}

/// Metadata attached to a registered modified-base symbol.
#[derive(Debug, Clone)]
pub struct ModSymbol {
    /// The single-character code used in sequence text (e.g. 'P' for pseudouridine).
    pub symbol: char,
    /// Canonical parent base this modification derives from (its digital code, `0..K`).
    pub parent: Dsq,
    /// Human-readable / MODOMICS short name (e.g. "pseudouridine", "m1A").
    pub name: String,
}

/// Specification for registering a modified base into an extended alphabet.
#[derive(Debug, Clone)]
pub struct ModSpec {
    /// Single-character sequence-text code for the modified base.
    pub symbol: char,
    /// Parent canonical base, as a character (one of A/C/G/U).
    pub parent: char,
    /// Short name (MODOMICS).
    pub name: String,
}

/// A biosequence alphabet (`ESL_ALPHABET`), possibly extended with modified bases.
#[derive(Debug, Clone)]
pub struct Alphabet {
    pub atype: AlphabetType,
    /// Number of canonical residues (4 for nucleic).
    pub k: usize,
    /// Total number of codes (canonical + gap + degenerate + modified + N/*/~).
    pub kp: usize,
    /// Symbol characters, indexed by digital code `[0..kp]`.
    pub sym: Vec<u8>,
    /// ASCII char -> digital code; [`ILLEGAL`] for unmapped chars.
    pub inmap: [Dsq; 128],
    /// `degen[code][0..k]`: which canonical residues each code includes.
    pub degen: Vec<Vec<bool>>,
    /// Number of canonical residues each code marginalizes over.
    pub ndegen: Vec<usize>,
    /// Complement code for each code (nucleic only).
    pub complement: Option<Vec<Dsq>>,
    /// Per-code modified-base metadata; `None` for ordinary codes.
    pub modinfo: Vec<Option<ModSymbol>>,
}

impl Alphabet {
    /// Standard RNA alphabet, faithful to Easel `create_rna()` (K=4, Kp=18).
    pub fn rna() -> Self {
        Self::nucleic(AlphabetType::Rna, &[])
    }

    /// Standard DNA alphabet, faithful to Easel `create_dna()` (K=4, Kp=18).
    pub fn dna() -> Self {
        Self::nucleic(AlphabetType::Dna, &[])
    }

    /// RNA alphabet **extended** with modified-base symbols. Each modified symbol is
    /// inserted into the degenerate region with a degeneracy vector pointing at its
    /// parent canonical base, so homology scoring is unchanged while the modified
    /// identity is preserved for annotation.
    pub fn rna_with_modifications(mods: &[ModSpec]) -> Result<Self, AlphabetError> {
        Self::nucleic_checked(AlphabetType::Rna, mods)
    }

    fn nucleic(atype: AlphabetType, mods: &[ModSpec]) -> Self {
        Self::nucleic_checked(atype, mods).expect("standard nucleic alphabet construction")
    }

    fn nucleic_checked(atype: AlphabetType, mods: &[ModSpec]) -> Result<Self, AlphabetError> {
        // Canonical residues and their IUPAC degeneracies, mirroring Easel exactly.
        // (RNA uses U; DNA uses T.)
        let (canon, complement_canon): (&[u8], &[u8]) = match atype {
            AlphabetType::Rna => (b"ACGU", b"UGCA"),
            AlphabetType::Dna => (b"ACGT", b"TGCA"),
            AlphabetType::Custom => {
                return Err(AlphabetError::Alphabet(
                    "custom nucleic not supported".into(),
                ))
            }
        };
        let k = canon.len(); // 4

        // (degenerate symbol, the canonical residues it covers) — Easel order RYMKSWHBVD.
        let u = if atype == AlphabetType::Rna {
            b'U'
        } else {
            b'T'
        };
        let degens: [(u8, Vec<u8>); 10] = [
            (b'R', vec![b'A', b'G']),
            (b'Y', vec![b'C', u]),
            (b'M', vec![b'A', b'C']),
            (b'K', vec![b'G', u]),
            (b'S', vec![b'C', b'G']),
            (b'W', vec![b'A', u]),
            (b'H', vec![b'A', b'C', u]),
            (b'B', vec![b'C', b'G', u]),
            (b'V', vec![b'A', b'C', b'G']),
            (b'D', vec![b'A', b'G', u]),
        ];

        // Build the symbol vector in canonical Easel order, splicing modified symbols
        // in just before the trailing N / * / ~ specials.
        let mut sym: Vec<u8> = Vec::new();
        sym.extend_from_slice(canon); // 0..K
        sym.push(b'-'); // K: gap
        for (d, _) in &degens {
            sym.push(*d); // K+1 .. K+10
        }
        let mod_start = sym.len();
        for m in mods {
            if !m.symbol.is_ascii() {
                return Err(AlphabetError::Alphabet(format!(
                    "modified symbol {:?} must be ASCII",
                    m.symbol
                )));
            }
            sym.push(m.symbol as u8);
        }
        sym.push(b'N'); // Kp-3: any
        sym.push(b'*'); // Kp-2: nonresidue
        sym.push(b'~'); // Kp-1: missing
        let kp = sym.len();

        // Canonical-residue index lookup for resolving parents/degeneracies.
        let canon_index: HashMap<u8, usize> =
            canon.iter().enumerate().map(|(i, &c)| (c, i)).collect();

        // Degeneracy matrix.
        let mut degen = vec![vec![false; k]; kp];
        let mut ndegen = vec![0usize; kp];
        for (i, _) in canon.iter().enumerate() {
            degen[i][i] = true;
            ndegen[i] = 1;
        }
        for (di, (_, residues)) in degens.iter().enumerate() {
            let code = k + 1 + di;
            for r in residues {
                let ci = canon_index[r];
                degen[code][ci] = true;
                ndegen[code] += 1;
            }
        }
        // "any" character (N) at Kp-3 covers all canonical residues.
        degen[kp - 3].fill(true);
        ndegen[kp - 3] = k;

        // Modified symbols: degeneracy = indicator of parent canonical base.
        let mut modinfo: Vec<Option<ModSymbol>> = vec![None; kp];
        for (mi, m) in mods.iter().enumerate() {
            let code = mod_start + mi;
            let parent_uc = (m.parent.to_ascii_uppercase()) as u8;
            // Map U/T to the canonical residue for this alphabet kind.
            let parent_byte = match (atype, parent_uc) {
                (AlphabetType::Rna, b'T') => b'U',
                (AlphabetType::Dna, b'U') => b'T',
                (_, b) => b,
            };
            let parent_idx = *canon_index.get(&parent_byte).ok_or_else(|| {
                AlphabetError::Alphabet(format!(
                    "modified symbol {:?} has non-canonical parent {:?}",
                    m.symbol, m.parent
                ))
            })?;
            degen[code][parent_idx] = true;
            ndegen[code] = 1;
            modinfo[code] = Some(ModSymbol {
                symbol: m.symbol,
                parent: parent_idx as Dsq,
                name: m.name.clone(),
            });
        }

        // Input map: base symbols, then Easel's synonyms, then case-insensitivity,
        // then modified symbols last so their chars win over any synonym collision.
        let mut inmap = [ILLEGAL; 128];
        for (x, &c) in sym.iter().enumerate() {
            inmap[c as usize] = x as Dsq;
        }
        // Synonyms (esl_alphabet_SetEquiv), matching create_rna/create_dna.
        let set_equiv = |from: u8, to: u8, inmap: &mut [Dsq; 128]| {
            if inmap[to as usize] != ILLEGAL {
                inmap[from as usize] = inmap[to as usize];
            }
        };
        match atype {
            AlphabetType::Rna => {
                set_equiv(b'T', b'U', &mut inmap);
                set_equiv(b'I', b'A', &mut inmap); // inosine ~ A, unless overridden below
            }
            AlphabetType::Dna => {
                set_equiv(b'U', b'T', &mut inmap);
                set_equiv(b'I', b'A', &mut inmap);
            }
            AlphabetType::Custom => {}
        }
        set_equiv(b'X', b'N', &mut inmap);
        set_equiv(b'_', b'-', &mut inmap);
        set_equiv(b'.', b'-', &mut inmap);
        // Case-insensitive (esl_alphabet_SetCaseInsensitive): lower follows upper.
        for c in b'A'..=b'Z' {
            let lc = c.to_ascii_lowercase();
            if inmap[c as usize] != ILLEGAL {
                inmap[lc as usize] = inmap[c as usize];
            }
        }
        // Modified symbols override synonyms/case for their exact characters.
        for (mi, m) in mods.iter().enumerate() {
            inmap[m.symbol as usize] = (mod_start + mi) as Dsq;
        }

        // Complement (set_complementarity), extended for modified symbols (-> parent's
        // complement, a canonical code) and any inserted positions.
        let mut complement = vec![0u8; kp];
        for (i, &cc) in complement_canon.iter().enumerate() {
            complement[i] = canon_index[&cc] as Dsq;
        }
        complement[k] = k as Dsq; // gap -> gap
                                  // Degenerate complements (R<->Y, M<->K, S/W self, H<->D, B<->V) by residue set.
        for di in 0..degens.len() {
            let code = k + 1 + di;
            // complement of a degenerate code = the code whose residue set is the
            // complement of this one. Compute by complementing residues then matching.
            let mut comp_residues: Vec<usize> = (0..k)
                .filter(|&x| degen[code][x])
                .map(|x| complement[x] as usize)
                .collect();
            comp_residues.sort_unstable();
            let mut found = code; // default to self
            for dj in 0..degens.len() {
                let other = k + 1 + dj;
                let mut other_res: Vec<usize> = (0..k).filter(|&x| degen[other][x]).collect();
                other_res.sort_unstable();
                if other_res == comp_residues {
                    found = other;
                    break;
                }
            }
            complement[code] = found as Dsq;
        }
        for (mi, _) in mods.iter().enumerate() {
            let code = mod_start + mi;
            let parent = modinfo[code].as_ref().unwrap().parent as usize;
            complement[code] = complement[parent]; // -> canonical complement of parent
        }
        complement[kp - 3] = (kp - 3) as Dsq; // N -> N
        complement[kp - 2] = (kp - 2) as Dsq; // * -> *
        complement[kp - 1] = (kp - 1) as Dsq; // ~ -> ~

        Ok(Alphabet {
            atype,
            k,
            kp,
            sym,
            inmap,
            degen,
            ndegen,
            complement: Some(complement),
            modinfo,
        })
    }

    /// Digital code for a canonical residue character (A/C/G/U), if it is one.
    pub fn canonical_code(&self, c: char) -> Option<Dsq> {
        if !c.is_ascii() {
            return None;
        }
        let code = self.inmap[c as usize];
        if (code as usize) < self.k {
            Some(code)
        } else {
            None
        }
    }

    /// Is `x` a canonical residue (`x < K`)?
    #[inline]
    pub fn is_canonical(&self, x: Dsq) -> bool {
        (x as usize) < self.k
    }

    /// Is `x` a residue (canonical, degenerate, modified, or N) — i.e. not gap/`*`/`~`?
    #[inline]
    pub fn is_residue(&self, x: Dsq) -> bool {
        let xi = x as usize;
        xi < self.k || (xi > self.k && xi < self.kp - 2)
    }

    /// Is `x` the gap code?
    #[inline]
    pub fn is_gap(&self, x: Dsq) -> bool {
        x as usize == self.k
    }

    /// Is `x` a registered modified base?
    #[inline]
    pub fn is_modified(&self, x: Dsq) -> bool {
        (x as usize) < self.kp && self.modinfo[x as usize].is_some()
    }

    /// Modified-base metadata for a code, if any.
    pub fn modinfo(&self, x: Dsq) -> Option<&ModSymbol> {
        self.modinfo.get(x as usize).and_then(|m| m.as_ref())
    }

    /// Digitize a text sequence into a 1-based digital sequence with sentinels at
    /// positions 0 and L+1 (`esl_abc_CreateDsq`). Residues live at `dsq[1..=L]`.
    pub fn digitize(&self, text: &str) -> Result<Vec<Dsq>, AlphabetError> {
        let bytes = text.as_bytes();
        let mut dsq = Vec::with_capacity(bytes.len() + 2);
        dsq.push(SENTINEL);
        for &b in bytes {
            if b.is_ascii_whitespace() {
                continue;
            }
            if b >= 128 {
                return Err(AlphabetError::Alphabet(format!(
                    "non-ASCII byte {b} in sequence"
                )));
            }
            let code = self.inmap[b as usize];
            if code == ILLEGAL {
                return Err(AlphabetError::Alphabet(format!(
                    "illegal symbol {:?} in sequence",
                    b as char
                )));
            }
            dsq.push(code);
        }
        dsq.push(SENTINEL);
        Ok(dsq)
    }

    /// Convert a 1-based digital sequence (with sentinels) back to text
    /// (`esl_abc_Textize`).
    pub fn textize(&self, dsq: &[Dsq]) -> String {
        dsq.iter()
            .filter(|&&x| x != SENTINEL)
            .map(|&x| self.sym[x as usize] as char)
            .collect()
    }

    /// Average score of code `x` over the canonical residues it marginalizes over
    /// (`esl_abc_FAvgScore`). For a canonical or modified symbol this returns the
    /// single residue's score; for a degenerate symbol it averages.
    #[inline]
    pub fn favg_score(&self, x: Dsq, sc: &[f32]) -> f32 {
        if !self.is_residue(x) {
            return 0.0;
        }
        let xi = x as usize;
        let mut result = 0.0f32;
        for (i, &s) in sc.iter().enumerate().take(self.k) {
            if self.degen[xi][i] {
                result += s;
            }
        }
        result / self.ndegen[xi] as f32
    }

    /// Reverse-complement a 1-based digital sequence in place (`esl_abc_revcomp`).
    pub fn revcomp(&self, dsq: &mut [Dsq]) -> Result<(), AlphabetError> {
        let comp = self
            .complement
            .as_ref()
            .ok_or_else(|| AlphabetError::Alphabet("alphabet has no complement".into()))?;
        // residues are at 1..=L; sentinels at 0 and L+1.
        if dsq.len() < 2 {
            return Ok(());
        }
        let l = dsq.len() - 2;
        // complement
        for x in dsq.iter_mut().skip(1).take(l) {
            *x = comp[*x as usize];
        }
        // reverse the residue span
        dsq[1..=l].reverse();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rna_layout_matches_easel() {
        let a = Alphabet::rna();
        assert_eq!(a.k, 4);
        assert_eq!(a.kp, 18);
        assert_eq!(
            a.sym.iter().map(|&b| b as char).collect::<String>(),
            "ACGU-RYMKSWHBVDN*~"
        );
        // canonical
        assert_eq!(a.inmap[b'A' as usize], 0);
        assert_eq!(a.inmap[b'U' as usize], 3);
        // T reads as U; lower-case follows upper.
        assert_eq!(a.inmap[b'T' as usize], 3);
        assert_eq!(a.inmap[b'a' as usize], 0);
        // degeneracy of R = {A,G}
        let r = a.inmap[b'R' as usize] as usize;
        assert_eq!(a.ndegen[r], 2);
        assert!(a.degen[r][0] && a.degen[r][2]);
    }

    #[test]
    fn digitize_roundtrip_with_sentinels() {
        let a = Alphabet::rna();
        let dsq = a.digitize("ACGURYN").unwrap();
        assert_eq!(dsq[0], SENTINEL);
        assert_eq!(*dsq.last().unwrap(), SENTINEL);
        assert_eq!(dsq.len(), 9);
        assert_eq!(a.textize(&dsq), "ACGURYN");
    }

    #[test]
    fn favg_marginalizes_degenerate() {
        let a = Alphabet::rna();
        let sc = [1.0f32, 2.0, 4.0, 8.0]; // A C G U
                                          // R = {A,G} -> (1+4)/2 = 2.5
        let r = a.inmap[b'R' as usize];
        assert!((a.favg_score(r, &sc) - 2.5).abs() < 1e-6);
        // canonical G -> 4
        assert!((a.favg_score(2, &sc) - 4.0).abs() < 1e-6);
    }

    #[test]
    fn revcomp_is_correct() {
        let a = Alphabet::rna();
        let mut dsq = a.digitize("AACGU").unwrap();
        a.revcomp(&mut dsq).unwrap();
        assert_eq!(a.textize(&dsq), "ACGUU");
    }

    #[test]
    fn modified_symbol_scores_like_parent() {
        // Pseudouridine 'P' (parent U), dihydrouridine 'D' would collide with degen D,
        // so use 'O' for a U-derived mod and 'L' for an A-derived mod in this test.
        let mods = vec![
            ModSpec {
                symbol: 'P',
                parent: 'U',
                name: "pseudouridine".into(),
            },
            ModSpec {
                symbol: 'L',
                parent: 'A',
                name: "m1A".into(),
            },
        ];
        let a = Alphabet::rna_with_modifications(&mods).unwrap();
        // Kp grew by 2.
        assert_eq!(a.kp, 20);
        assert_eq!(a.k, 4);

        let sc = [1.0f32, 2.0, 4.0, 8.0];
        let p = a.inmap[b'P' as usize];
        let l = a.inmap[b'L' as usize];
        assert!(a.is_modified(p));
        assert!(a.is_residue(p));
        // P scores exactly like U (8.0); L scores exactly like A (1.0).
        assert!((a.favg_score(p, &sc) - 8.0).abs() < 1e-6);
        assert!((a.favg_score(l, &sc) - 1.0).abs() < 1e-6);
        assert_eq!(a.modinfo(p).unwrap().name, "pseudouridine");
        assert_eq!(a.modinfo(p).unwrap().parent, 3); // U

        // Digitize a sequence containing a modified base; it round-trips with identity.
        let dsq = a.digitize("ACGUP").unwrap();
        assert_eq!(a.textize(&dsq), "ACGUP");
    }

    #[test]
    fn extended_alphabet_preserves_canonical_scoring() {
        // The key invariant: adding modified symbols does not perturb canonical/degenerate
        // scoring relative to the stock alphabet.
        let base = Alphabet::rna();
        let mods = vec![ModSpec {
            symbol: 'P',
            parent: 'U',
            name: "Y".into(),
        }];
        let ext = Alphabet::rna_with_modifications(&mods).unwrap();
        let sc = [1.5f32, -0.3, 2.2, 0.7];
        for c in [b'A', b'C', b'G', b'U', b'R', b'Y', b'N'] {
            let xb = base.inmap[c as usize];
            let xe = ext.inmap[c as usize];
            assert!(
                (base.favg_score(xb, &sc) - ext.favg_score(xe, &sc)).abs() < 1e-6,
                "score mismatch for {}",
                c as char
            );
        }
    }
}
