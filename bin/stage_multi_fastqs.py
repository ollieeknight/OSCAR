#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

from stage_fastqs import stage_modality

MODALITY_FEATURE = {
    "gex":    "Gene Expression",
    "adt":    "Antibody Capture",
    "hto":    "Antibody Capture",
    "vdj_t":  "VDJ-T",
    "vdj_b":  "VDJ-B",
    "crispr": "CRISPR Guide Capture",
}


def stage(library_id, min_reads):
    with open("config_header.txt") as fh:
        config_header = fh.read().strip()

    with open("flex_samples.txt") as fh:
        flex_samples_content = fh.read().strip()

    work_dir  = str(Path.cwd())
    lib_lines = []

    for mod, feat in MODALITY_FEATURE.items():
        staged = sorted(
            p for p in Path(".").glob(f"fastqs/{mod}/run_*/*")
            if p.name != "NO_FILE"
        )
        if not staged:
            continue

        try:
            lane, sample_id = stage_modality(staged, Path("fastq_all") / mod, min_reads)
        except ValueError as exc:
            print(f"[cellranger_multi] ERROR: {mod}: {exc}", file=sys.stderr)
            sys.exit(1)

        if lane > 0 and sample_id:
            lib_lines.append(f"{sample_id},{work_dir}/fastq_all/{mod},{feat}")

    cfg = config_header
    if lib_lines:
        cfg += "\n\n[libraries]\nfastq_id,fastqs,feature_types\n"
        cfg += "\n".join(lib_lines)
    if flex_samples_content:
        cfg += "\n" + flex_samples_content
    cfg = re.sub(r'\n{3,}', '\n\n', cfg)

    with open("multi_config.csv", "w") as fh:
        fh.write(cfg + "\n")

    print("[cellranger_multi] Config written:", file=sys.stderr)
    with open("multi_config.csv") as fh:
        print(fh.read(), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library-id", required=True,
                    help="library id embedded in the multi config")
    ap.add_argument("--min-reads", type=int, required=True,
                    help="lanes with fewer reads than this are skipped")
    args = ap.parse_args()
    stage(args.library_id, args.min_reads)


if __name__ == "__main__":
    main()
