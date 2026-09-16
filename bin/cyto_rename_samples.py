#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path


def rename(samples_file, cyto_out):
    samples_path = Path(samples_file)
    cyto_out     = Path(cyto_out)
    counts_dir   = cyto_out / "counts"

    if samples_path.name == "NO_FILE" or not counts_dir.exists():
        print("[cyto_rename] no samples file or no counts dir — skipping rename")
    else:
        mapping = {}
        with open(samples_path) as f:
            for row in csv.DictReader(f):
                sample = row['sample_id'].strip()
                for bc in row['probe_barcode_ids'].strip().replace('|', ',').split(','):
                    mapping[bc.strip()] = sample

        for bc_dir in sorted(counts_dir.iterdir()):
            if bc_dir.name in mapping:
                target = bc_dir.parent / mapping[bc_dir.name]
                shutil.move(str(bc_dir), str(target))
                print(f"[cyto_rename] {bc_dir.name} → {mapping[bc_dir.name]}")
            else:
                print(f"[cyto_rename] no mapping for {bc_dir.name} — left as-is")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples-file", required=True,
                    help="CSV mapping sample_id to probe_barcode_ids, or NO_FILE")
    ap.add_argument("--cyto-out", required=True,
                    help="cyto output directory containing counts/")
    args = ap.parse_args()
    rename(args.samples_file, args.cyto_out)


if __name__ == "__main__":
    main()
