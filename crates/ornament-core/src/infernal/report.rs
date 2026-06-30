//! Tab-separated hit report for the native scanner.
//!
//! A wider, still tab-separated alternative to the bare coordinate dump: it surfaces each hit's
//! model coordinate span, GC, the model accession, the per-hit significance marker, and the
//! target description. It deliberately omits Infernal `--tblout`'s always-constant columns
//! (`mdl`, `trunc`, `pass`, `bias`) and its fixed-width padding — the values are the same, the
//! layout is clean TSV.

use super::CMHit;

/// E-value at or below which a hit is marked significant (`!`), mirroring cmsearch's default
/// inclusion threshold (`--incE 0.01`); above it the marker is `?`.
const INC_THRESHOLD: f64 = 0.01;

/// Tab-separated column header for [`write_tsv`].
pub const TSV_HEADER: &str = "target_name\ttarget_start\ttarget_end\tstrand\tquery_name\tquery_accession\tmdl_from\tmdl_to\tgc\tscore\te_value\tinc\tdescription";

/// Render hits as a tab-separated table with a header row (trailing newline included).
pub fn write_tsv(hits: &[CMHit]) -> String {
    let mut out = String::with_capacity(TSV_HEADER.len() + hits.len() * 64);
    out.push_str(TSV_HEADER);
    out.push('\n');
    write_tsv_rows(&mut out, hits);
    out
}

/// Append the tab-separated data rows (no header) for `hits` to `out`. Used by [`write_tsv`] and by
/// the streaming scan path, which writes the header once then appends each batch as it completes.
pub fn write_tsv_rows(out: &mut String, hits: &[CMHit]) {
    for h in hits {
        // Uncalibrated models report a NaN E-value; show `-` and leave the hit unmarked.
        let inc = if h.e_value <= INC_THRESHOLD { '!' } else { '?' };
        let evalue = if h.e_value.is_nan() {
            "-".to_string()
        } else {
            format!("{:.2e}", h.e_value)
        };
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.1}\t{}\t{}\t{}\n",
            h.target_name,
            h.target_start,
            h.target_end,
            h.strand,
            h.query_name,
            h.query_accession.as_deref().unwrap_or("-"),
            h.mdl_from,
            h.mdl_to,
            h.gc_content,
            h.score,
            evalue,
            inc,
            h.description.as_deref().unwrap_or("-"),
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(strand: char, e_value: f64) -> CMHit {
        CMHit {
            target_name: "chr1".to_string(),
            target_start: 41,
            target_end: 111,
            strand,
            query_name: "tRNA".to_string(),
            score: 87.2,
            e_value,
            gc_content: 0.58,
            mdl_from: 1,
            mdl_to: 71,
            query_accession: Some("RF00005".to_string()),
            description: None,
        }
    }

    #[test]
    fn header_and_row_layout() {
        let out = write_tsv(&[hit('+', 5e-30)]);
        let mut lines = out.lines();
        assert_eq!(lines.next().unwrap(), TSV_HEADER);
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(row.len(), 13, "one field per header column");
        assert_eq!(row[0], "chr1");
        assert_eq!(row[3], "+");
        assert_eq!(row[5], "RF00005");
        assert_eq!(row[6], "1"); // mdl_from
        assert_eq!(row[7], "71"); // mdl_to
        assert_eq!(row[8], "0.58"); // gc, %.2
        assert_eq!(row[9], "87.2"); // score, %.1
        assert_eq!(row[10], "5.00e-30"); // e_value, %.2e
        assert_eq!(row[11], "!"); // significant (<= 0.01)
        assert_eq!(row[12], "-"); // absent description
    }

    #[test]
    fn significance_marker_and_nan_evalue() {
        let out = write_tsv(&[hit('-', 0.5), hit('+', f64::NAN)]);
        let rows: Vec<&str> = out.lines().skip(1).collect();
        // E-value above the inclusion threshold -> `?`.
        assert_eq!(rows[0].split('\t').nth(11).unwrap(), "?");
        // Uncalibrated (NaN) -> `-` E-value, unmarked.
        let nan_row: Vec<&str> = rows[1].split('\t').collect();
        assert_eq!(nan_row[10], "-");
        assert_eq!(nan_row[11], "?");
    }
}
