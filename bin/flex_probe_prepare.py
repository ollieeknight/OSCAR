#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def prepare(standard_probe_csv, custom_probe_csv):
    standard_rows = []
    std_path = Path(standard_probe_csv)
    if std_path.name != "NO_FILE":
        with open(std_path) as f:
            lines = [l for l in f if not l.startswith('#')]
            standard_rows = list(csv.DictReader(lines))

    custom_rows = []
    cust_path = Path(custom_probe_csv)
    if cust_path.name != "NO_FILE":
        with open(cust_path) as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                custom_rows.append({
                    'gene_id':   row['gene_id'],
                    'probe_seq': row['probe_seq'],
                    'probe_id':  row['probe_id'],
                    'included':  'TRUE' if row['included'].upper() == 'TRUE' else 'FALSE',
                    'region':    row.get('region', 'spliced'),
                    'gene_name': row.get('gene_name', row['gene_id']),
                })

    seen = {}
    for r in standard_rows + custom_rows:
        seen[r['probe_id']] = r
    merged = list(seen.values())

    cr_fields = ['gene_id', 'probe_seq', 'probe_id', 'included']
    if merged and 'region' in merged[0]:
        cr_fields.append('region')
    with open("merged_probes_cellranger.csv", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cr_fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged)

    with open("merged_probes_cyto.tsv", "w") as f:
        for r in merged:
            if r.get('included', 'TRUE').upper() == 'TRUE':
                gene_name = r.get('gene_name', r['gene_id'])
                f.write(f"{r['probe_id']}\t{gene_name}\t{r['probe_seq']}\n")

    print(f"Merged {len(merged)} probes ({len(standard_rows)} standard + {len(custom_rows)} custom)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--standard-probe-csv", required=True,
                    help="standard 10x probe CSV, or NO_FILE")
    ap.add_argument("--custom-probe-csv", required=True,
                    help="custom OSCAR-format probe CSV, or NO_FILE")
    args = ap.parse_args()
    prepare(args.standard_probe_csv, args.custom_probe_csv)


if __name__ == "__main__":
    main()
