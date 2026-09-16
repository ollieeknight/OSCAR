#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

from stage_fastqs import count_reads, pair_reads


def stage(min_reads):
    staged = sorted(p for p in Path(".").glob("fastqs/gex/run_*/*") if p.name != "NO_FILE")

    read_pairs, unmatched = pair_reads(staged)
    if unmatched:
        print(f"[cyto_flex] ERROR: unmatched reads: {unmatched}", file=sys.stderr)
        sys.exit(1)

    pairs = []
    for r1, r2 in read_pairs:
        n = count_reads(r1, min_reads)
        if n >= min_reads:
            pairs.extend([str(r1), str(r2)])
        else:
            print(f"[cyto_flex] skip {r1.parent.name}: {n} reads (<{min_reads})", file=sys.stderr)

    if not pairs:
        print("[cyto_flex] ERROR: no FASTQ pairs passed read threshold", file=sys.stderr)
        sys.exit(1)

    with open("fastq_pairs.txt", "w") as f:
        f.write(" ".join(pairs))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-reads", type=int, required=True,
                    help="pairs with fewer reads than this are skipped")
    args = ap.parse_args()
    stage(args.min_reads)


if __name__ == "__main__":
    main()
