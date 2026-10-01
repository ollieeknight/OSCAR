#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# Column names are bcl-convert's Reports/Demultiplex_Stats.csv and Top_Unknown_Barcodes.csv.
def read_demultiplex_stats(path):
    per_lane = defaultdict(dict)
    for row in _rows(path):
        lane, sample = row["Lane"], row["SampleID"]
        per_lane[lane][sample] = per_lane[lane].get(sample, 0) + int(row["# Reads"])
    return dict(per_lane)


def read_top_unknown(path):
    # index2 is absent on single-index runs.
    return [(row["Lane"], row["index"], row.get("index2") or "", int(row["# Reads"]))
            for row in _rows(path)]


def load_kit_indexes(indexes_dir):
    table = {}
    for csv_path in sorted(Path(indexes_dir).glob("*.csv")):
        with open(csv_path, newline="", encoding="utf-8-sig") as fh:
            for parts in csv.reader(fh):
                if not parts or not parts[0].strip():
                    continue
                name = parts[0].strip()
                if name.lower().startswith("index_name"):
                    continue
                for cell in parts[1:]:
                    seq = cell.strip().upper()
                    if seq and set(seq) <= set("ACGTN"):
                        table.setdefault(seq, name)
    return table


def _median(values):
    s = sorted(values)
    if not s:
        return 0
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2


def summarise(per_lane):
    rows = []
    for lane in sorted(per_lane, key=lambda x: (len(x), x)):
        counts = per_lane[lane]
        total = sum(counts.values())
        for sample, reads in sorted(counts.items(), key=lambda kv: -kv[1]):
            pct = (100.0 * reads / total) if total else 0.0
            rows.append(
                {"lane": lane, "library": sample, "reads": reads, "percent_of_lane": f"{pct:.4f}"}
            )
    return rows


def flowcell_overview(per_lane):
    totals = defaultdict(int)
    lanes_by_sample = defaultdict(list)
    for lane in sorted(per_lane, key=lambda x: (len(x), x)):
        for sample, reads in per_lane[lane].items():
            totals[sample] += reads
            lanes_by_sample[sample].append(lane)
    grand_total = sum(totals.values())
    rows = []
    for sample, reads in sorted(totals.items(), key=lambda kv: -kv[1]):
        pct = (100.0 * reads / grand_total) if grand_total else 0.0
        rows.append(
            {
                "name": sample,
                "read_number": reads,
                "percent_of_flowcell": f"{pct:.4f}",
                "lanes": ",".join(lanes_by_sample[sample]),
            }
        )
    return rows


def find_dropouts(per_lane, ratio):
    warnings = []
    for lane in sorted(per_lane, key=lambda x: (len(x), x)):
        declared = {s: n for s, n in per_lane[lane].items() if s.lower() != "undetermined"}
        if len(declared) < 2:
            continue
        med = _median(list(declared.values()))
        if med <= 0:
            continue
        for sample, reads in sorted(declared.items(), key=lambda kv: kv[1]):
            if reads * ratio < med:
                warnings.append(
                    {
                        "lane": lane,
                        "type": "dropout",
                        "library": sample,
                        "index": "",
                        "reads": reads,
                        "percent_of_lane": "",
                        "detail": (
                            f"declared library has {reads} reads against a lane median of {med:.0f}"
                            f" -- likely wrong index in metadata.csv or failed library"
                        ),
                    }
                )
    return warnings


def find_undeclared(unknown, per_lane, kit_indexes, min_pct):
    warnings = []
    for lane, i7, i5, reads in unknown:
        total = sum(per_lane.get(lane, {}).values())
        if not total:
            continue
        pct = 100.0 * reads / total
        if pct < min_pct:
            continue
        hits = [kit_indexes[s] for s in (i7.upper(), i5.upper()) if s in kit_indexes]
        if not hits:
            continue
        label = "+".join(filter(None, [i7, i5]))
        warnings.append(
            {
                "lane": lane,
                "type": "undeclared",
                "library": "",
                "index": label,
                "reads": reads,
                "percent_of_lane": f"{pct:.4f}",
                "detail": (
                    f"unknown barcode at {pct:.2f}% of lane reads matches kit index "
                    f"{'/'.join(sorted(set(hits)))} -- library sequenced but absent from metadata.csv"
                ),
            }
        )
    return sorted(warnings, key=lambda w: -w["reads"])


def _write(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True, help="Demultiplex_Stats.csv")
    ap.add_argument("--unknown", required=True, help="Top_Unknown_Barcodes.csv")
    ap.add_argument("--indexes", required=True, help="assets/indexes directory")
    ap.add_argument("--prefix", required=True, help="output file prefix (run name)")
    ap.add_argument("--dropout-ratio", type=float, required=True)
    ap.add_argument("--unknown-pct", type=float, required=True)
    args = ap.parse_args()

    per_lane = read_demultiplex_stats(args.stats)
    summary = summarise(per_lane)
    _write(
        f"{args.prefix}_demux_summary.csv",
        ["lane", "library", "reads", "percent_of_lane"],
        summary,
    )

    _write(
        f"{args.prefix}_flow_cell_overview.csv",
        ["name", "read_number", "percent_of_flowcell", "lanes"],
        flowcell_overview(per_lane),
    )

    warnings = find_dropouts(per_lane, args.dropout_ratio) + find_undeclared(
        read_top_unknown(args.unknown),
        per_lane,
        load_kit_indexes(args.indexes),
        args.unknown_pct,
    )
    _write(
        f"{args.prefix}_demux_warnings.csv",
        ["lane", "type", "library", "index", "reads", "percent_of_lane", "detail"],
        warnings,
    )

    for w in warnings:
        print(
            f"WARNING: [{w['type']}] lane {w['lane']} "
            f"{w['library'] or w['index']}: {w['detail']}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
