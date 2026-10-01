#!/usr/bin/env python3
"""Stage FASTQs for cellranger multi, cellranger-atac or cyto.

Nextflow stages inputs as fastqs/<modality>/run_NNN/*. Read sets with fewer than
--min-reads reads are dropped; the rest are renamed to <sample>_S1_L<n>_<read>_001.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

MULTI_FEATURES = {
    "gex":    "Gene Expression",
    "adt":    "Antibody Capture",
    "hto":    "Antibody Capture",
    "vdj_t":  "VDJ-T",
    "vdj_b":  "VDJ-B",
    "crispr": "CRISPR Guide Capture",
}


def flowcell_of(fastq):
    with gzip.open(fastq, "rt") as fh:
        header = fh.readline().strip()
    fields = header.lstrip("@").split(":")
    if len(fields) < 4:
        sys.exit(f"ERROR: unreadable FASTQ header in {fastq}: {header!r}")
    return tuple(fields[:3])


def count_reads(fastq, limit):
    try:
        with gzip.open(fastq, "rt") as fh:
            return sum(1 for i, _ in zip(range(limit * 4), fh) if i % 4 == 0)
    except (OSError, EOFError):
        return 0


def read_sets(modality, reads, min_reads):
    """Group a modality's FASTQs by flowcell and lane, keeping sets with enough reads."""
    sets = {}
    for p in Path("fastqs", modality).glob("run_*/*"):
        m = re.search(r"_(R[123])_", p.name)
        if m and m.group(1) in reads:
            key = (flowcell_of(p), re.sub(r"_R[123]_", "_R_", p.name))
            sets.setdefault(key, {})[m.group(1)] = p

    kept = []
    for key in sorted(sets):
        missing = [r for r in reads if r not in sets[key]]
        if missing:
            sys.exit(f"ERROR: {key[1]} has no {'/'.join(missing)}")
        n = count_reads(sets[key][reads[0]], min_reads)
        if n < min_reads:
            print(f"skip {key[1]} from flowcell {key[0][2]}: {n} reads (<{min_reads})", file=sys.stderr)
        else:
            kept.append(sets[key])
    return kept


def rename(sets, out_dir):
    out_dir.mkdir(parents=True)
    sample = re.sub(r"_S\d+_.*", "", next(iter(sets[0].values())).name)
    for lane, reads in enumerate(sets, 1):
        for read, p in reads.items():
            p.rename(out_dir / f"{sample}_S1_L{lane:03d}_{read}_001.fastq.gz")
    return sample


def multi(min_reads):
    libraries = []
    for modality, feature in MULTI_FEATURES.items():
        sets = read_sets(modality, ["R1", "R2"], min_reads)
        if sets:
            out_dir = Path.cwd() / "staged" / modality
            libraries.append(f"{rename(sets, out_dir)},{out_dir},{feature}")

    config = Path("config_header.txt").read_text().strip()
    if libraries:
        config += "\n\n[libraries]\nfastq_id,fastqs,feature_types\n" + "\n".join(libraries)
    samples = Path("flex_samples.txt").read_text().strip()
    if samples:
        config += "\n\n" + samples
    Path("multi_config.csv").write_text(config + "\n")
    print(config, file=sys.stderr)


def atac(min_reads):
    sets = read_sets("atac", ["R1", "R2", "R3"], min_reads)
    if not sets:
        sys.exit("ERROR: no ATAC read sets passed the read threshold")
    rename(sets, Path("staged/atac"))


def cyto(min_reads):
    sets = read_sets("gex", ["R1", "R2"], min_reads)
    if not sets:
        sys.exit("ERROR: no GEX read pairs passed the read threshold")
    Path("fastq_pairs.txt").write_text(" ".join(f"{s['R1']} {s['R2']}" for s in sets))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", choices=["multi", "atac", "cyto"])
    ap.add_argument("--min-reads", type=int, default=10000)
    args = ap.parse_args()
    {"multi": multi, "atac": atac, "cyto": cyto}[args.target](args.min_reads)


if __name__ == "__main__":
    main()
