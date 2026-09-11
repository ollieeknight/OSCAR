#!/usr/bin/env python3
"""Summarise a bcl-convert Reports/ directory and flag suspect demultiplexing.

Imported by tests; invoked by DEMUX_QC via PATH.

Reads bcl-convert's own Demultiplex_Stats.csv and Top_Unknown_Barcodes.csv --
the stats it already writes for free -- instead of re-deriving read counts by
decompressing multi-GB Undetermined FASTQ files.

Two outputs:
  <run>_demux_summary.csv   library, reads, percent of lane total (Undetermined
                            included, so percentages sum to 100).
  <run>_demux_warnings.csv  one row per suspect finding, empty (header only)
                            when the demux looks clean.

Checks:
  dropout   a declared library whose read count is far below the median of its
            siblings in the same lane. A library that failed to demux shows up
            orders of magnitude down, not merely low, hence a ratio threshold.
  undeclared  an index in Top_Unknown_Barcodes.csv above a share of lane reads
            that matches a known kit index. An unknown barcode that is also a
            real kit index means a library was sequenced but left out of
            metadata.csv -- the case that has to be caught, as opposed to
            ordinary free-floating sequencing noise.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# A declared library this many times below its lane's median is treated as
# failed rather than merely under-loaded. Matches the ratio that surfaced the
# real dropouts when these runs were diagnosed by hand.
DEFAULT_DROPOUT_RATIO = 100.0
# An unknown barcode below this share of lane reads is noise, not a library.
DEFAULT_UNKNOWN_PCT = 1.0


def _rows(path):
    """Read a bcl-convert CSV, tolerating a UTF-8 BOM and blank trailing lines."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return [r for r in csv.DictReader(fh) if any(v.strip() for v in r.values() if v)]


def _int(value):
    """bcl-convert writes read counts plain, but tolerate thousands separators."""
    return int(str(value).strip().replace(",", "") or 0)


def _col(row, *names):
    """Fetch the first present column. Header spelling drifts between versions."""
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return None


def read_demultiplex_stats(path):
    """-> {lane: {sample_id: reads}}. Lanes are kept separate; a library absent
    from one lane's pool is normal and must not look like a dropout there."""
    per_lane = defaultdict(dict)
    for row in _rows(path):
        sample = _col(row, "SampleID", "Sample_ID", "Sample ID")
        if sample is None:
            continue
        lane = str(_col(row, "Lane") or "1")
        reads = _int(_col(row, "# Reads", "Reads", "NumberOfReads") or 0)
        # Same sample can appear once per index row; accumulate rather than overwrite.
        per_lane[lane][sample] = per_lane[lane].get(sample, 0) + reads
    return dict(per_lane)


def read_top_unknown(path):
    """-> [(lane, index, index2, reads)] from Top_Unknown_Barcodes.csv."""
    out = []
    for row in _rows(path):
        lane = str(_col(row, "Lane") or "1")
        i7 = (_col(row, "index", "Index", "Barcode") or "").strip()
        i5 = (_col(row, "index2", "Index2") or "").strip()
        out.append((lane, i7, i5, _int(_col(row, "# Reads", "Reads", "Count") or 0)))
    return out


def load_kit_indexes(indexes_dir):
    """-> {sequence: kit_code} for every sequence in assets/indexes/*.csv.

    Kit CSVs are headerless or headered depending on the kit, and carry a
    varying number of sequence columns (4 for single-index, 3 for dual), so
    every non-name cell that looks like a barcode is indexed."""
    table = {}
    for csv_path in sorted(Path(indexes_dir).glob("*.csv")):
        with open(csv_path, newline="", encoding="utf-8-sig") as fh:
            for parts in csv.reader(fh):
                if not parts or not parts[0].strip():
                    continue
                name = parts[0].strip()
                if name.lower().startswith("index_name"):
                    continue  # header row, present in some kits only
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
    """-> summary rows, percent taken against the lane total including Undetermined."""
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


def find_dropouts(per_lane, ratio=DEFAULT_DROPOUT_RATIO):
    """Declared libraries sitting far below their lane's median sibling.

    The median is taken over the real libraries only -- Undetermined is often
    the largest or smallest entry in a lane and would skew the baseline."""
    warnings = []
    for lane in sorted(per_lane, key=lambda x: (len(x), x)):
        declared = {s: n for s, n in per_lane[lane].items() if s.lower() != "undetermined"}
        if len(declared) < 2:
            continue  # no siblings to compare against; nothing to conclude
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
                            f"declared library has {reads} reads, "
                            f"{med / reads:.0f}x below lane median {med:.0f} "
                            f"-- likely wrong index in metadata.csv or failed library"
                        )
                        if reads
                        else (
                            f"declared library has 0 reads "
                            f"(lane median {med:.0f}) "
                            f"-- likely wrong index in metadata.csv or failed library"
                        ),
                    }
                )
    return warnings


def find_undeclared(unknown, per_lane, kit_indexes, min_pct=DEFAULT_UNKNOWN_PCT):
    """Unknown barcodes above min_pct of lane reads that match a known kit index."""
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
            continue  # unknown but not a kit index: sequencing noise, not a library
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats", required=True, help="Demultiplex_Stats.csv")
    ap.add_argument("--unknown", help="Top_Unknown_Barcodes.csv")
    ap.add_argument("--indexes", help="assets/indexes directory")
    ap.add_argument("--prefix", required=True, help="output file prefix (run name)")
    ap.add_argument("--dropout-ratio", type=float, default=DEFAULT_DROPOUT_RATIO)
    ap.add_argument("--unknown-pct", type=float, default=DEFAULT_UNKNOWN_PCT)
    args = ap.parse_args(argv)

    per_lane = read_demultiplex_stats(args.stats)
    summary = summarise(per_lane)
    _write(
        f"{args.prefix}_demux_summary.csv",
        ["lane", "library", "reads", "percent_of_lane"],
        summary,
    )

    warnings = find_dropouts(per_lane, args.dropout_ratio)
    if args.unknown and args.indexes and Path(args.unknown).exists():
        warnings += find_undeclared(
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
    # Warnings are advisory: a flagged demux must stay visible without failing
    # the run, since the FASTQ are still valid output.
    return 0


if __name__ == "__main__":
    sys.exit(main())
