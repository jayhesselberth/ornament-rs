//! Parser for the INFERNAL1/a ASCII `.cm` covariance-model format.
//!
//! Faithful port of `read_asc_1p1_cm` in `src/cm_file.c`. The on-disk file stores
//! log-odds *scores*; we reconstruct probabilities with `ascii2prob` exactly:
//! `p = (s == "*") ? 0 : 2^score · null`. Log-odds scores (`tsc`/`esc`) are recomputed
//! later by [`crate::config`]. The embedded HMMER3 filter HMM is retained as raw text
//! for Phase 5; the CM E-value `ECM*` lines are parsed into [`ExpInfo`].

use std::path::Path;
use std::sync::Arc;

use easel_rs::Alphabet;

use crate::model::{
    derive_unique_state_code, exp_mode, nd, node_code, st, state_code, Cm, CmFlags, ExpInfo,
    NodeAnnot, MAXCONNECT,
};
use crate::InfernalError;

/// `ascii2prob`: a `.cm` score token -> probability. `2^score · null`, or 0 for `*`.
#[inline]
fn ascii2prob(s: &str, null: f32) -> Result<f32, InfernalError> {
    if s == "*" {
        return Ok(0.0);
    }
    let score: f64 = s
        .parse()
        .map_err(|_| InfernalError::Parse(format!("invalid score token {s:?}")))?;
    // p = 2^score · null. Since score/log2(e) = score·ln2, exp(score/log2(e)) = 2^score.
    // 1.4426950408889634 == log2(e), matching Infernal's `ascii2prob` constant.
    Ok((score / 1.4426950408889634).exp() as f32 * null)
}

fn opt_i32(tok: &str) -> Option<i32> {
    if tok == "-" {
        None
    } else {
        tok.parse().ok()
    }
}

fn opt_char(tok: &str) -> Option<char> {
    if tok == "-" {
        None
    } else {
        tok.chars().next()
    }
}

/// Parse a `.cm` file from disk.
pub fn parse_cm_file<P: AsRef<Path>>(path: P) -> Result<Cm, InfernalError> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|e| InfernalError::Io(format!("{}: {e}", path.as_ref().display())))?;
    parse_cm_str(&text)
}

/// Parse a `.cm` model from text.
pub fn parse_cm_str(text: &str) -> Result<Cm, InfernalError> {
    let mut lines = text.lines().peekable();

    // First line: format tag.
    let first = lines
        .next()
        .ok_or_else(|| InfernalError::Parse("empty CM file".into()))?;
    if !first.starts_with("INFERNAL1/a") {
        return Err(InfernalError::Parse(format!(
            "expected INFERNAL1/a tag, found {:?}",
            first.split_whitespace().next().unwrap_or("")
        )));
    }

    // --- Header section, until a bare "CM" line. ---
    let mut h = Header::default();
    for line in lines.by_ref() {
        let mut it = line.split_whitespace();
        let Some(tag) = it.next() else { continue };
        if tag == "CM" {
            break;
        }
        let rest: Vec<&str> = it.collect();
        h.apply(tag, &rest)?;
    }

    let abc = Arc::new(Alphabet::rna());
    let k = abc.k;
    let m = h.states.ok_or_else(|| InfernalError::Parse("missing STATES".into()))?;
    let nodes = h.nodes.ok_or_else(|| InfernalError::Parse("missing NODES".into()))?;
    let clen_hdr = h.clen.ok_or_else(|| InfernalError::Parse("missing CLEN".into()))?;
    let null = h
        .null
        .clone()
        .unwrap_or_else(|| vec![1.0 / k as f32; k]);

    // --- Body section: node lines and state lines, until "//". ---
    let mut cm = Cm {
        name: h.name.clone().unwrap_or_default(),
        acc: h.acc.clone(),
        desc: h.desc.clone(),
        rf: None,
        consensus: None,
        map: None,
        nseq: h.nseq.unwrap_or(0),
        eff_nseq: h.eff_nseq.unwrap_or(0.0),
        null,
        m,
        clen: 0, // recomputed below, then checked against header
        nodes,
        w: h.w.unwrap_or(0),
        sttype: vec![0; m],
        ndidx: vec![0; m],
        stid: vec![0; m],
        cfirst: vec![0; m],
        cnum: vec![0; m],
        plast: vec![0; m],
        pnum: vec![0; m],
        nodemap: vec![0; nodes],
        ndtype: vec![0; nodes],
        node_annot: vec![NodeAnnot::default(); nodes],
        dmin1: vec![0; m],
        dmin2: vec![0; m],
        dmax1: vec![0; m],
        dmax2: vec![0; m],
        t: vec![Vec::new(); m],
        e: vec![Vec::new(); m],
        tsc: vec![Vec::new(); m],
        esc: vec![Vec::new(); m],
        begin: vec![0.0; m],
        end: vec![0.0; m],
        pbegin: h.pbegin.unwrap_or(0.05),
        pend: h.pend.unwrap_or(0.05),
        el_selfsc: h.el_selfsc.unwrap_or(0.0),
        ga: h.ga,
        tc: h.tc,
        nc: h.nc,
        exp: h.exp,
        fp7_gfmu: h.fp7_gfmu,
        fp7_gflambda: h.fp7_gflambda,
        fp7_text: None,
        flags: h.flags,
        abc,
    };

    let mut cur_nd: usize = 0;
    let mut v: usize = 0;
    let mut hmm_text: Option<String> = None;

    while let Some(line) = lines.next() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.is_empty() {
            continue;
        }
        if toks[0] == "//" {
            // Anything after the CM's closing // is the embedded HMMER3 filter HMM.
            let mut buf = String::new();
            for l in lines.by_ref() {
                buf.push_str(l);
                buf.push('\n');
            }
            if !buf.trim().is_empty() {
                hmm_text = Some(buf);
            }
            break;
        }
        if toks[0] == "[" {
            // v is the index of the first state that will belong to this node.
            parse_node_line(&mut cm, &toks, &mut cur_nd, v)?;
        } else {
            parse_state_line(&mut cm, &toks, cur_nd, v, k)?;
            v += 1;
        }
    }

    cm.fp7_text = hmm_text;

    // --- Validate against header. ---
    if v != m {
        return Err(InfernalError::Parse(format!(
            "read {v} states, header declared STATES {m}"
        )));
    }
    if cm.clen != clen_hdr {
        return Err(InfernalError::Parse(format!(
            "computed consensus length {} != header CLEN {clen_hdr}",
            cm.clen
        )));
    }

    Ok(cm)
}

fn parse_node_line(
    cm: &mut Cm,
    toks: &[&str],
    cur_nd: &mut usize,
    first_state: usize,
) -> Result<(), InfernalError> {
    // [ <NDTYPE> <nd> ] <lmap> <rmap> <lcons> <rcons> <lrf> <rrf>
    if toks.len() < 4 {
        return Err(InfernalError::Parse(format!("short node line: {toks:?}")));
    }
    let ndtype = node_code(toks[1])
        .ok_or_else(|| InfernalError::Parse(format!("invalid node type {:?}", toks[1])))?;
    let nd_idx: usize = toks[2]
        .parse()
        .map_err(|_| InfernalError::Parse(format!("invalid node index {:?}", toks[2])))?;
    if nd_idx >= cm.nodes {
        return Err(InfernalError::Parse(format!(
            "node index {nd_idx} exceeds NODES {}",
            cm.nodes
        )));
    }
    cm.ndtype[nd_idx] = ndtype;
    match ndtype {
        nd::MATP => cm.clen += 2,
        nd::MATL | nd::MATR => cm.clen += 1,
        _ => {}
    }
    // nodemap[nd] = index of this node's first state (the next state line read).
    cm.nodemap[nd_idx] = first_state;
    *cur_nd = nd_idx;

    // Annotation columns after the closing ']' (index 3): map(2) cons(2) rf(2).
    let a = &toks[4..];
    let mut ann = NodeAnnot::default();
    if a.len() >= 2 {
        ann.lmap = opt_i32(a[0]);
        ann.rmap = opt_i32(a[1]);
    }
    if a.len() >= 4 {
        ann.lcons = opt_char(a[2]);
        ann.rcons = opt_char(a[3]);
    }
    if a.len() >= 6 {
        ann.lrf = opt_char(a[4]);
        ann.rrf = opt_char(a[5]);
    }
    cm.node_annot[nd_idx] = ann;
    Ok(())
}

fn parse_state_line(
    cm: &mut Cm,
    toks: &[&str],
    cur_nd: usize,
    v: usize,
    k: usize,
) -> Result<(), InfernalError> {
    // <sttype> <v> <plast> <pnum> <cfirst> <cnum> <dmin2> <dmin1> <dmax1> <dmax2> <t...> <e...>
    if v >= cm.m {
        return Err(InfernalError::Parse(format!(
            "more state lines than STATES ({})",
            cm.m
        )));
    }
    if toks.len() < 10 {
        return Err(InfernalError::Parse(format!("short state line: {toks:?}")));
    }
    let sttype = state_code(toks[0])
        .ok_or_else(|| InfernalError::Parse(format!("invalid state type {:?}", toks[0])))?;
    let declared_v: usize = toks[1]
        .parse()
        .map_err(|_| InfernalError::Parse(format!("invalid state index {:?}", toks[1])))?;
    if declared_v != v {
        return Err(InfernalError::Parse(format!(
            "state index {declared_v} out of order (expected {v})"
        )));
    }
    let parse_i = |t: &str, what: &str| -> Result<i32, InfernalError> {
        t.parse()
            .map_err(|_| InfernalError::Parse(format!("invalid {what} {t:?} on state {v}")))
    };

    cm.sttype[v] = sttype;
    cm.plast[v] = parse_i(toks[2], "plast")?;
    cm.pnum[v] = parse_i(toks[3], "pnum")?;
    cm.cfirst[v] = parse_i(toks[4], "cfirst")?;
    cm.cnum[v] = parse_i(toks[5], "cnum")?;
    cm.dmin2[v] = parse_i(toks[6], "dmin2")?;
    cm.dmin1[v] = parse_i(toks[7], "dmin1")?;
    cm.dmax1[v] = parse_i(toks[8], "dmax1")?;
    cm.dmax2[v] = parse_i(toks[9], "dmax2")?;

    let mut idx = 10;

    // Transitions: cnum of them, except B_st (which carries no transition probs;
    // its cfirst/cnum encode the two bifurcation children).
    if sttype != st::B {
        let cnum = cm.cnum[v].max(0) as usize;
        if cnum > MAXCONNECT {
            return Err(InfernalError::Parse(format!(
                "state {v} cnum {cnum} exceeds MAXCONNECT {MAXCONNECT}"
            )));
        }
        let mut tv = Vec::with_capacity(cnum);
        for x in 0..cnum {
            let tok = toks.get(idx + x).ok_or_else(|| {
                InfernalError::Parse(format!("too few transition fields on state {v}"))
            })?;
            tv.push(ascii2prob(tok, 1.0)?);
        }
        idx += cnum;
        cm.t[v] = tv;
    }

    // Emissions.
    if matches!(sttype, st::ML | st::MR | st::IL | st::IR) {
        let mut ev = Vec::with_capacity(k);
        for x in 0..k {
            let tok = toks
                .get(idx + x)
                .ok_or_else(|| InfernalError::Parse(format!("too few emissions on state {v}")))?;
            ev.push(ascii2prob(tok, cm.null[x])?);
        }
        cm.e[v] = ev;
    } else if sttype == st::MP {
        let mut ev = Vec::with_capacity(k * k);
        for x in 0..k {
            for y in 0..k {
                let tok = toks.get(idx + x * k + y).ok_or_else(|| {
                    InfernalError::Parse(format!("too few pair emissions on state {v}"))
                })?;
                ev.push(ascii2prob(tok, cm.null[x] * cm.null[y])?);
            }
        }
        cm.e[v] = ev;
    }

    cm.ndidx[v] = cur_nd;
    cm.stid[v] = derive_unique_state_code(cm.ndtype[cur_nd], sttype).ok_or_else(|| {
        InfernalError::Parse(format!(
            "illegal state type {sttype} in node type {}",
            cm.ndtype[cur_nd]
        ))
    })?;
    Ok(())
}

/// Accumulator for header key/value lines.
#[derive(Default)]
struct Header {
    name: Option<String>,
    acc: Option<String>,
    desc: Option<String>,
    states: Option<usize>,
    nodes: Option<usize>,
    clen: Option<usize>,
    w: Option<i32>,
    nseq: Option<i32>,
    eff_nseq: Option<f32>,
    null: Option<Vec<f32>>,
    pbegin: Option<f32>,
    pend: Option<f32>,
    el_selfsc: Option<f32>,
    ga: Option<f32>,
    tc: Option<f32>,
    nc: Option<f32>,
    fp7_gfmu: Option<f64>,
    fp7_gflambda: Option<f64>,
    exp: [ExpInfo; exp_mode::NMODES],
    flags: CmFlags,
}

impl Header {
    fn apply(&mut self, tag: &str, rest: &[&str]) -> Result<(), InfernalError> {
        let first = rest.first().copied();
        let yes = first == Some("yes");
        match tag {
            "NAME" => self.name = first.map(str::to_string),
            "ACC" => self.acc = first.map(str::to_string),
            "DESC" => self.desc = Some(rest.join(" ")),
            "STATES" => self.states = first.and_then(|s| s.parse().ok()),
            "NODES" => self.nodes = first.and_then(|s| s.parse().ok()),
            "CLEN" => self.clen = first.and_then(|s| s.parse().ok()),
            "W" => self.w = first.and_then(|s| s.parse().ok()),
            "RF" => self.flags.rf = yes,
            "CONS" => self.flags.cons = yes,
            "MAP" => self.flags.map = yes,
            "NSEQ" => self.nseq = first.and_then(|s| s.parse().ok()),
            "EFFN" => self.eff_nseq = first.and_then(|s| s.parse().ok()),
            "PBEGIN" => self.pbegin = first.and_then(|s| s.parse().ok()),
            "PEND" => self.pend = first.and_then(|s| s.parse().ok()),
            "ELSELF" => self.el_selfsc = first.and_then(|s| s.parse().ok()),
            "GA" => {
                self.ga = first.and_then(|s| s.parse().ok());
                self.flags.ga = self.ga.is_some();
            }
            "TC" => {
                self.tc = first.and_then(|s| s.parse().ok());
                self.flags.tc = self.tc.is_some();
            }
            "NC" => {
                self.nc = first.and_then(|s| s.parse().ok());
                self.flags.nc = self.nc.is_some();
            }
            "NULL" => {
                // NULL stores log-odds vs uniform; reconstruct probabilities.
                let k = rest.len();
                let mut nv = Vec::with_capacity(k);
                for s in rest {
                    nv.push(ascii2prob(s, 1.0 / k as f32)?);
                }
                self.null = Some(nv);
            }
            "EFP7GF" => {
                self.fp7_gfmu = rest.first().and_then(|s| s.parse().ok());
                self.fp7_gflambda = rest.get(1).and_then(|s| s.parse().ok());
            }
            _ if tag.starts_with("ECM") => {
                let mode = match tag {
                    "ECMLC" => exp_mode::CM_LC,
                    "ECMGC" => exp_mode::CM_GC,
                    "ECMLI" => exp_mode::CM_LI,
                    "ECMGI" => exp_mode::CM_GI,
                    _ => return Err(InfernalError::Parse(format!("unknown ECM tag {tag}"))),
                };
                if rest.len() < 6 {
                    return Err(InfernalError::Parse(format!("short {tag} line")));
                }
                let p = |i: usize| -> f64 { rest[i].parse().unwrap_or(0.0) };
                self.exp[mode] = ExpInfo {
                    lambda: p(0),
                    mu_extrap: p(1),
                    mu_orig: p(2),
                    dbsize: p(3),
                    nrandhits: rest[4].parse().unwrap_or(0),
                    tailp: p(5),
                    is_valid: true,
                };
            }
            _ => {} // ignore COM/DATE/WBETA/QDBBETA*/N2OMEGA/N3OMEGA/CKSUM/etc.
        }
        Ok(())
    }
}
