#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path

from stage_fastqs import count_reads, flowcell_of


def stage(sample_id, min_reads):
    staged = sorted(
        p for p in Path(".").glob("fastqs/atac/run_*/*")
        if p.name != "NO_FILE"
    )

    runs = {}
    for p in staged:
        m = re.search(r'_(R[123])_', p.name)
        if not m:
            continue                   # index reads (I1/I2) are not inputs
        key = (flowcell_of(p), re.sub(r'_R[123]_', '_R_', p.name))
        runs.setdefault(key, {})[m.group(1)] = p

    if not runs:
        print("[cellranger_atac] ERROR: no R1/R2/R3 FASTQs staged", file=sys.stderr)
        sys.exit(1)

    final_dir = Path("fastq_all/atac")
    final_dir.mkdir(parents=True, exist_ok=True)
    lane = 0

    for key in sorted(runs):
        reads = runs[key]
        fq_group = reads["R1"].name if "R1" in reads else str(key)
        missing = [t for t in ("R1", "R2", "R3") if t not in reads]
        if missing:
            print(f"[cellranger_atac] ERROR: {fq_group} missing {missing}", file=sys.stderr)
            sys.exit(1)
        n = count_reads(reads["R1"], min_reads)
        if n < min_reads:
            print(f"[cellranger_atac] skip {fq_group}: {n} reads", file=sys.stderr)
            continue
        lane += 1
        for tag, p in sorted(reads.items()):
            p.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_{tag}_001.fastq.gz")

    if lane == 0:
        print("[cellranger_atac] ERROR: no valid FASTQ lanes found", file=sys.stderr)
        sys.exit(1)

    with open("fastq_dir.txt", "w") as fh:
        fh.write(str(final_dir.resolve()))

    return lane


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True,
                    help="sample id to embed in the renamed FASTQ filenames")
    ap.add_argument("--min-reads", type=int, required=True,
                    help="lanes with fewer reads than this are skipped")
    args = ap.parse_args()
    stage(args.sample_id, args.min_reads)


if __name__ == "__main__":
    main()
