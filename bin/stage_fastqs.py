#!/usr/bin/env python3
"""Pair staged FASTQs by flowcell and rename them into per-modality lanes.

Imported by tests; invoked by CELLRANGER_MULTI via `python3 -c` / PATH.

Nextflow's `stageAs: "run_???/*"` gives every input file its own run_??? dir,
so the dir name records input-list position, not which flowcell a file came
from. Two flowcells routinely emit identically-named FASTQs, so neither the
filename nor the dir can pair R1 with R2 -- only the read header can.
"""
import gzip
import re
import subprocess
import sys
from pathlib import Path


def flowcell_of(fastq, _open=gzip.open):
    """(instrument, run, flowcell) parsed from the first read header."""
    with _open(fastq, "rt") as fh:
        header = fh.readline().strip()
    fields = header.lstrip("@").split(":")
    if len(fields) < 4:
        raise ValueError(f"unreadable FASTQ header in {fastq}: {header!r}")
    return tuple(fields[0:3])


def pair_reads(staged, header_of=flowcell_of):
    """Pair each R1 with the R2 from the same flowcell.

    Returns (pairs, unmatched). Pairs are ordered by original filename then
    flowcell, so lane numbering is stable across resumes.
    """
    def key(p):
        return (header_of(p), re.sub(r"_R[12]_", "_R_", p.name))

    r1 = {key(p): p for p in staged if re.search(r"_R1_", p.name)}
    r2 = {key(p): p for p in staged if re.search(r"_R2_", p.name)}

    unmatched = [r1[k].name for k in sorted(r1) if k not in r2]
    unmatched += [r2[k].name for k in sorted(r2) if k not in r1]

    pairs = [(r1[k], r2[k]) for k in sorted(r1, key=lambda k: (k[1], k[0]))
             if k in r2]
    return pairs, unmatched


def count_reads(r1, min_reads):
    """Reads in the first `min_reads` records; cheap placeholder detection."""
    result = subprocess.run(
        f"zcat {r1} | head -n {min_reads * 4} | awk 'NR%4==1' | wc -l",
        shell=True, capture_output=True, text=True,
    )
    return int(result.stdout.strip() or "0")


def stage_modality(staged, final_dir, min_reads, counter=count_reads):
    """Rename paired FASTQs into final_dir as {sample}_S1_L00N_R[12]_001.

    Returns (lanes_written, sample_id). Placeholder pairs below min_reads are
    skipped without consuming a lane number.
    """
    pairs, unmatched = pair_reads(staged)
    if unmatched:
        raise ValueError(f"unmatched reads: {unmatched}")
    if not pairs:
        return 0, None

    final_dir.mkdir(parents=True, exist_ok=True)
    lane, sample_id = 0, None
    for r1, r2 in pairs:
        n = counter(r1, min_reads)
        if n < min_reads:
            print(f"[stage_fastqs] skip {r1.parent.name}: {n} reads", file=sys.stderr)
            continue
        lane += 1
        sample_id = re.sub(r"_S\d+.*", "", r1.name)
        r1.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R1_001.fastq.gz")
        r2.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R2_001.fastq.gz")
    return lane, sample_id
