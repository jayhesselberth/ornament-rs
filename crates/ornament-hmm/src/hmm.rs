//! Parser for the HMMER3/f profile-HMM format (the embedded CM filter HMM, `fp7`).
//!
//! The on-disk format stores parameters as `-ln(prob)` (`*` = probability 0). We convert
//! back to probabilities; [`crate::profile`] turns those into the log-odds search profile.
//! Field order mirrors HMMER's `p7_hmmfile`: per node a match-emission line, an
//! insert-emission line, and a 7-value transition line
//! (`m→m m→i m→d i→m i→i d→m d→d`), preceded by a COMPO line and the node-0 (begin) block.

use crate::HmmerError;

/// Transition indices into the 7-value per-node transition array.
pub mod tr {
    pub const MM: usize = 0;
    pub const MI: usize = 1;
    pub const MD: usize = 2;
    pub const IM: usize = 3;
    pub const II: usize = 4;
    pub const DM: usize = 5;
    pub const DD: usize = 6;
}

/// Gumbel/exponential calibration for a filter stage (`mu`, `lambda`).
#[derive(Debug, Clone, Copy, Default)]
pub struct EvParam {
    pub mu: f64,
    pub lambda: f64,
}

/// A Plan7 profile HMM (`P7_HMM`), parameters in probability space.
#[derive(Debug, Clone)]
pub struct P7Hmm {
    pub name: String,
    pub acc: Option<String>,
    /// Model length = number of match states.
    pub m: usize,
    /// Alphabet size K (4 for RNA).
    pub k: usize,
    pub maxl: i32,
    /// Match emission probs `mat[k][a]`, k in 1..=M (mat[0] unused / COMPO-adjacent).
    pub mat: Vec<Vec<f32>>,
    /// Insert emission probs `ins[k][a]`, k in 0..=M.
    pub ins: Vec<Vec<f32>>,
    /// Transition probs `t[k][0..7]`, k in 0..=M (node 0 = begin transitions).
    pub t: Vec<[f32; 7]>,
    /// Overall composition (COMPO line), if present.
    pub compo: Option<Vec<f32>>,
    pub msv: EvParam,
    pub vit: EvParam,
    pub fwd: EvParam,
}

/// `-ln(prob)` token → probability. `*` → 0.
fn p(tok: &str) -> Result<f32, HmmerError> {
    if tok == "*" {
        return Ok(0.0);
    }
    let v: f64 = tok
        .parse()
        .map_err(|_| HmmerError::Parse(format!("invalid prob token {tok:?}")))?;
    Ok((-v).exp() as f32)
}

/// Parse a HMMER3/f HMM from text (the section between `HMMER3/f` and the closing `//`).
pub fn parse_p7_hmm(text: &str) -> Result<P7Hmm, HmmerError> {
    let mut lines = text.lines().peekable();
    let first = lines
        .next()
        .ok_or_else(|| HmmerError::Parse("empty HMM".into()))?;
    if !first.starts_with("HMMER3/") {
        return Err(HmmerError::Parse(format!(
            "expected HMMER3/ tag, found {first:?}"
        )));
    }

    let mut name = String::new();
    let mut acc = None;
    let mut m = 0usize;
    let mut k = 4usize;
    let mut maxl = 0i32;
    let (mut msv, mut vit, mut fwd) = (EvParam::default(), EvParam::default(), EvParam::default());

    // Header until the "HMM" line (which is followed by the alphabet + transition headers).
    for line in lines.by_ref() {
        let mut it = line.split_whitespace();
        let Some(tag) = it.next() else { continue };
        match tag {
            "HMM" => break,
            "NAME" => name = it.next().unwrap_or("").to_string(),
            "ACC" => acc = it.next().map(str::to_string),
            "LENG" => m = it.next().and_then(|s| s.parse().ok()).unwrap_or(0),
            "MAXL" => maxl = it.next().and_then(|s| s.parse().ok()).unwrap_or(0),
            "ALPH" => {
                let a = it.next().unwrap_or("");
                k = if a.eq_ignore_ascii_case("amino") {
                    20
                } else {
                    4
                };
            }
            "STATS" => {
                // STATS LOCAL <MSV|VITERBI|FORWARD> <mu> <lambda>
                let _local = it.next();
                let which = it.next().unwrap_or("");
                let mu = it.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let lambda = it.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let ev = EvParam { mu, lambda };
                match which {
                    "MSV" => msv = ev,
                    "VITERBI" => vit = ev,
                    "FORWARD" => fwd = ev,
                    _ => {}
                }
            }
            _ => {}
        }
    }
    if m == 0 {
        return Err(HmmerError::Parse("missing LENG".into()));
    }

    // The "HMM" line already carries the alphabet header (`HMM  A C G U`); the next line is
    // the transition header (`m->m m->i ...`), which we skip. COMPO/node-0 follow.
    lines.next(); // transition header

    // COMPO line (optional), then node-0 block (insert-0 emissions, node-0 transitions).
    let mut compo = None;
    let mut mat: Vec<Vec<f32>> = vec![vec![0.0; k]; m + 1];
    let mut ins: Vec<Vec<f32>> = vec![vec![0.0; k]; m + 1];
    let mut t: Vec<[f32; 7]> = vec![[0.0; 7]; m + 1];

    let read_k = |line: &str, n: usize| -> Result<Vec<f32>, HmmerError> {
        let toks: Vec<&str> = line.split_whitespace().collect();
        toks.iter().take(n).map(|s| p(s)).collect()
    };
    let read_t = |line: &str| -> Result<[f32; 7], HmmerError> {
        let toks: Vec<&str> = line.split_whitespace().collect();
        let mut out = [0.0f32; 7];
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = p(toks.get(i).copied().unwrap_or("*"))?;
        }
        Ok(out)
    };

    let l0 = lines
        .next()
        .ok_or_else(|| HmmerError::Parse("missing COMPO/node0".into()))?;
    if l0.trim_start().starts_with("COMPO") {
        compo = Some(read_k(l0.trim_start().trim_start_matches("COMPO"), k)?);
        // node-0 insert emissions, then node-0 transitions
        ins[0] = read_k(
            lines
                .next()
                .ok_or_else(|| HmmerError::Parse("missing ins0".into()))?,
            k,
        )?;
        t[0] = read_t(
            lines
                .next()
                .ok_or_else(|| HmmerError::Parse("missing t0".into()))?,
        )?;
    } else {
        // No COMPO: l0 is the node-0 insert emissions.
        ins[0] = read_k(l0, k)?;
        t[0] = read_t(
            lines
                .next()
                .ok_or_else(|| HmmerError::Parse("missing t0".into()))?,
        )?;
    }

    // Nodes 1..=M: match-emission line (k values + annotations), insert line, transition line.
    for node in 1..=m {
        let mline = lines
            .next()
            .ok_or_else(|| HmmerError::Parse(format!("missing match line for node {node}")))?;
        // First token is the node index; the next k tokens are emissions.
        let mut it = mline.split_whitespace();
        let idx: usize = it
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| HmmerError::Parse(format!("bad node index at {node}")))?;
        if idx != node {
            return Err(HmmerError::Parse(format!(
                "node index {idx} out of order (expected {node})"
            )));
        }
        let mut me = Vec::with_capacity(k);
        for _ in 0..k {
            me.push(p(it.next().ok_or_else(|| {
                HmmerError::Parse(format!("too few emissions at node {node}"))
            })?)?);
        }
        mat[node] = me;
        ins[node] = read_k(
            lines
                .next()
                .ok_or_else(|| HmmerError::Parse("missing ins".into()))?,
            k,
        )?;
        t[node] = read_t(
            lines
                .next()
                .ok_or_else(|| HmmerError::Parse("missing t".into()))?,
        )?;
    }

    Ok(P7Hmm {
        name,
        acc,
        m,
        k,
        maxl,
        mat,
        ins,
        t,
        compo,
        msv,
        vit,
        fwd,
    })
}
