#!/usr/bin/env python3
"""Turn `cellranger multi` per-sample metrics into a MultiQC custom-content table.

Imported by tests; invoked by CELLRANGER_MQC via PATH.

MultiQC's built-in `cellranger` module only parses web summaries from
`cellranger count` and `cellranger vdj` (it matches on
'"subcommand":"count"' / '"subcommand":"vdj"'). This pipeline runs
`cellranger multi`, whose web_summary.html declares subcommand "multi" and is
therefore never picked up. Rather than leave the metrics unreported, the
per-sample metrics_summary.csv files are pivoted into a custom-content table.

Input layout (cellranger multi):
    <library>/outs/per_sample_outs/<sample>/metrics_summary.csv

That CSV is long format, one metric per row:
    Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value

Only the headline per-sample Gene Expression metrics are surfaced; the full
file stays on disk for anyone who needs the rest. A wide table with every
metric of every library type is unreadable in a report.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# Headline metrics, in display order. Matched case-insensitively against
# "Metric Name". Kept short deliberately: this table is for spotting an outlier
# library at a glance, not for replacing metrics_summary.csv.
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
    """'1,234' -> 1234.0, '92.5%' -> 92.5, 'N/A' -> None.

    Values are kept numeric where possible so MultiQC sorts and colours the
    column rather than treating it as free text."""
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
    """-> {metric_name: value} for the per-sample Gene Expression metrics.

    A metric name can repeat across library types and grouping levels, so rows
    are filtered to the sample-level Gene Expression view — the same numbers the
    web summary shows under 'Cells'. Where a metric appears only at library
    level, the library-level value is used as a fallback so the column is not
    silently empty."""
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
    """-> [(sample_id, path)] for every per_sample_outs metrics_summary.csv.

    The sample id is the per_sample_outs subdirectory name, which is what
    cellranger uses to identify a demultiplexed sample."""
    out = []
    for path in sorted(Path(root).rglob("per_sample_outs/*/metrics_summary.csv")):
        out.append((path.parent.name, path))
    return out


def build_table(entries):
    """-> (ordered metric names actually present, {sample: {metric: value}})."""
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
    """Write a MultiQC custom-content CSV with a header comment block."""
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
        return 0  # nothing to report is not an error; MultiQC just omits the section

    present, rows = build_table(entries)
    if not rows:
        print("No recognised metrics found in metrics_summary.csv", file=sys.stderr)
        return 0
    write_mqc(args.out, present, rows)
    print(f"Wrote {args.out}: {len(rows)} sample(s), {len(present)} metric(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
