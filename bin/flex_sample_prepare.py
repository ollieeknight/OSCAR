#!/usr/bin/env python3

import argparse
import csv


def prepare(probe_barcodes_ref, samples_file):
    bc_map = {}
    with open(probe_barcodes_ref) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            variant, canonical, bc_id = parts[0], parts[1], parts[2]
            bc_map.setdefault(bc_id, []).append((variant, canonical))

    used_bcs = set()
    with open(samples_file) as f:
        for row in csv.DictReader(f):
            for bc in row['probe_barcode_ids'].strip().replace('|', ',').split(','):
                used_bcs.add(bc.strip())

    with open("cyto_probe_barcodes.txt", "w") as out:
        for bc_id in sorted(used_bcs):
            entries = bc_map.get(bc_id, [])
            if not entries:
                print(f"WARNING: no sequences found for {bc_id} in probe barcode reference")
            for variant, canonical in entries:
                out.write(f"{variant}\t{canonical}\t{bc_id}\n")

    print(f"Wrote barcode sequences for {len(used_bcs)} probe barcodes")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--probe-barcodes-ref", required=True,
                    help="10x probe barcode reference file")
    ap.add_argument("--samples-file", required=True,
                    help="CSV with a probe_barcode_ids column")
    args = ap.parse_args()
    prepare(args.probe_barcodes_ref, args.samples_file)


if __name__ == "__main__":
    main()
