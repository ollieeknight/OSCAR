#!/usr/bin/env python3
"""Write a bcl-convert v2 SampleSheet.csv for one demux group.

--specs is a TSV of sample_id, i7, i5 (- if none), mask for 4-read runs, mask for 3-read runs.
A mask segment ending in * takes whatever cycles the read has left.
"""

import argparse
import csv
import re
import sys
import xml.etree.ElementTree as ET


def expand(mask, cycles):
    parts = mask.split(";")
    if len(parts) > len(cycles):
        sys.exit(f"ERROR: mask {mask} has {len(parts)} segments but the run has {len(cycles)} reads")
    out = []
    for part, n in zip(parts, cycles):
        if not part.endswith("*"):
            out.append(part)
            continue
        base = part[:-1]
        rest = n - sum(int(x) for x in re.findall(r"\d+", base))
        if rest < 0:
            sys.exit(f"ERROR: mask segment {part} needs more than the {n} cycles in its read")
        out.append(f"{base}{rest}" if rest else base[:-1])
    return ";".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-info", required=True)
    ap.add_argument("--specs", required=True)
    ap.add_argument("--dual", action="store_true")
    args = ap.parse_args()

    cycles = [int(r.get("NumCycles")) for r in ET.parse(args.run_info).iter("Read")]
    with open(args.specs, newline="") as fh:
        specs = [row for row in csv.reader(fh, delimiter="\t") if row]
    rows = [(sid, i7, i5.strip("-"), expand(oc4 if len(cycles) == 4 else oc3, cycles))
            for sid, i7, i5, oc4, oc3 in specs]
    # bcl-convert rejects OverrideCycles given both globally and per sample.
    per_sample = len({mask for *_, mask in rows}) > 1

    names = ["Read1Cycles", "Index1Cycles"] + (["Index2Cycles"] if len(cycles) == 4 else []) + ["Read2Cycles"]
    lines = ["[Header]", "FileFormatVersion,2", "", "[Reads]"]
    lines += [f"{name},{n}" for name, n in zip(names, cycles)]
    lines += ["", "[BCLConvert_Settings]"]
    if not per_sample:
        lines.append(f"OverrideCycles,{rows[0][3]}")
    lines.append("BarcodeMismatchesIndex1,1")
    if args.dual and not per_sample:
        lines.append("BarcodeMismatchesIndex2,1")
    lines += ["", "[BCLConvert_Data]"]
    lines.append(",".join(["Sample_ID", "Index"] + (["Index2"] if args.dual else []) + (["OverrideCycles"] if per_sample else [])))
    for sid, i7, i5, mask in rows:
        lines.append(",".join([sid, i7] + ([i5] if args.dual else []) + ([mask] if per_sample else [])))

    with open("SampleSheet.csv", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
