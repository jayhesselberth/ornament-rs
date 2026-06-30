#!/usr/bin/env python3
"""Build the Sprinzl-position reference covariance model.

tRNAscan-SE's `*_sprinzl.cm` reference models were never publicly released
(confirmed by the authors in UCSC-LoweLab/tRNAscan-SE#14), so we build our own from
the clover project's Sprinzl ground truth. The result is an all-match covariance
model whose consensus column *i* is exactly Sprinzl label *i* (variable positions
absent in a given tRNA align as deletions).

Inputs:
  - clover global-coords TSV (per-residue trna_id, seq_index, sprinzl_label,
    global_index, residue), e.g.
    clover/inst/extdata/sprinzl/sacCer_global_coords.tsv.gz

Outputs (written next to this repo's test fixtures):
  - sprinzl_euk.cm            (cmbuild --hand, clen == #global columns)
  - sprinzl_gcol_label.tsv    (consensus/global column -> Sprinzl label)

Requires Infernal's `cmbuild` on PATH (use: pixi run -e dev python scripts/build-sprinzl-cm.py).
"""
import csv
import collections
import gzip
import os
import subprocess
import sys

# Canonical Sprinzl base pairs (from tRNAscan-SE SprinzlPos.pm).
SPRINZL_PAIRS = [
    ("1", "72"), ("2", "71"), ("3", "70"), ("4", "69"), ("5", "68"), ("6", "67"), ("7", "66"),
    ("10", "25"), ("11", "24"), ("12", "23"), ("13", "22"),
    ("27", "43"), ("28", "42"), ("29", "41"), ("30", "40"), ("31", "39"),
    ("e11", "e21"), ("e12", "e22"), ("e13", "e23"), ("e14", "e24"),
    ("e15", "e25"), ("e16", "e26"), ("e17", "e27"),
    ("49", "65"), ("50", "64"), ("51", "63"), ("52", "62"), ("53", "61"),
]


def main(coords_tsv: str, out_dir: str):
    open_fn = gzip.open if coords_tsv.endswith(".gz") else open
    rows = list(csv.DictReader(open_fn(coords_tsv, "rt"), delimiter="\t"))

    by_trna = collections.defaultdict(dict)
    gcol_label = {}
    for r in rows:
        g = int(r["global_index"])
        res = {"T": "U"}.get(r["residue"].upper(), r["residue"].upper())
        res = res if res in "ACGU" else "N"
        by_trna[r["trna_id"]][g] = res
        gcol_label[g] = r["sprinzl_label"]
    ncols = max(gcol_label)

    label_gcol = {v: k for k, v in gcol_label.items()}
    ss = ["."] * ncols
    for a, b in SPRINZL_PAIRS:
        if a in label_gcol and b in label_gcol:
            ss[label_gcol[a] - 1] = "<"
            ss[label_gcol[b] - 1] = ">"

    seed = os.path.join(out_dir, "sprinzl_seed.sto")
    with open(seed, "w") as out:
        out.write("# STOCKHOLM 1.0\n\n")
        for t, gd in by_trna.items():
            out.write(f"{t:40s} " + "".join(gd.get(g, '-') for g in range(1, ncols + 1)) + "\n")
        out.write(f"{'#=GC SS_cons':40s} " + "".join(ss) + "\n")
        out.write(f"{'#=GC RF':40s} " + ("x" * ncols) + "\n")  # all-match
        out.write("//\n")

    cm = os.path.join(out_dir, "sprinzl_euk.cm")
    subprocess.run(["cmbuild", "--hand", "-F", cm, seed], check=True)

    with open(os.path.join(out_dir, "sprinzl_gcol_label.tsv"), "w") as out:
        out.write("global_col\tsprinzl_label\n")
        for g in range(1, ncols + 1):
            out.write(f"{g}\t{gcol_label.get(g, '?')}\n")
    print(f"built {cm} (clen={ncols}) from {len(by_trna)} tRNAs")


if __name__ == "__main__":
    coords = sys.argv[1] if len(sys.argv) > 1 else (
        "/beevol/home/jhessel/devel/rnabioco/clover/inst/extdata/sprinzl/"
        "sacCer_global_coords.tsv.gz"
    )
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(__file__), "..", "crates", "ornament-scfg", "tests", "data"
    )
    main(coords, out)
