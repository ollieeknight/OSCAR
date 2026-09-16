#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

WANTED = [
    "Cells",
    "Mean reads per cell",
    "Median genes per cell",
    "Median UMI counts per cell",
    "Confidently mapped reads in cells",
    "Fraction reads in cells",
    "Estimated number of cells",
    "Number of reads",
    "Valid barcodes",
    "Sequencing saturation",
]


def _num(value):
    v = str(value).strip().replace(",", "")
    if not v or v.upper() in ("N/A", "NA", "NONE"):
        return None
    pct = v.endswith("%")
    if pct:
        v = v[:-1]
    try:
        return float(v)
    except ValueError:
        return None


def parse_metrics(path):
    primary, fallback = {}, {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("Metric Name") or "").strip()
            if not name:
                continue
            value = row.get("Metric Value")
            lib = (row.get("Library Type") or "").strip().lower()
            cat = (row.get("Category") or "").strip().lower()
            if lib and lib != "gene expression":
                continue
            if cat == "cells":
                primary.setdefault(name, value)
            else:
                fallback.setdefault(name, value)
    merged = dict(fallback)
    merged.update(primary)
    return merged


def find_metric_files(root):
    out = []
    for path in sorted(Path(root).rglob("per_sample_outs/*/metrics_summary.csv")):
        out.append((path.parent.name, path))
    return out


def build_table(entries):
    rows = {}
    for sample, path in entries:
        metrics = parse_metrics(path)
        picked = {}
        lowered = {k.lower(): (k, v) for k, v in metrics.items()}
        for want in WANTED:
            hit = lowered.get(want.lower())
            if hit:
                picked[want] = hit[1]
        if picked:
            rows[sample] = picked
    present = [m for m in WANTED if any(m in r for r in rows.values())]
    return present, rows


def write_mqc(path, present, rows, section_id="oscar_cellranger"):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# id: {section_id}\n")
        fh.write("# section_name: 'Cell Ranger multi'\n")
        fh.write(
            "# description: 'Per-sample metrics from cellranger multi "
            "(per_sample_outs/*/metrics_summary.csv). MultiQC's built-in "
            "cellranger module only reads count/vdj web summaries, so these "
            "are parsed directly.'\n"
        )
        fh.write("# plot_type: 'table'\n")
        fh.write("# pconfig:\n")
        fh.write("#     id: 'oscar_cellranger_table'\n")
        w = csv.writer(fh)
        w.writerow(["Sample"] + present)
        for sample in sorted(rows):
            w.writerow([sample] + [rows[sample].get(m, "") for m in present])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=".", help="directory to search for cellranger outs")
    ap.add_argument("--out", required=True, help="output *_mqc.csv path")
    args = ap.parse_args(argv)

    entries = find_metric_files(args.root)
    if not entries:
        print("No per_sample_outs/*/metrics_summary.csv found", file=sys.stderr)
        return 0

    present, rows = build_table(entries)
    if not rows:
        print("No recognised metrics found in metrics_summary.csv", file=sys.stderr)
        return 0
    write_mqc(args.out, present, rows)
    print(f"Wrote {args.out}: {len(rows)} sample(s), {len(present)} metric(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
