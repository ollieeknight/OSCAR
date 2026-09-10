#!/usr/bin/env python3
"""Self-checks for bin/stage_fastqs.py. Run: python3 tests/test_stage_fastqs.py

Regression target: cellranger TXRNGR10013, seen on R502 when four flowcells
were merged. fastq_all/hto/..._L003_R1 came from one flowcell and _L003_R2
from another, so cellranger compared headers 222KF7NNX vs 222KG5FNX at line 4
and aborted. Half the lanes paired correctly, which is why it surfaced late.
"""
import gzip
import io
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
from stage_fastqs import flowcell_of, pair_reads, stage_modality  # noqa: E402


def write_fastq(path, instrument, run, flowcell, read):
    """One-record gzipped FASTQ with a realistic Illumina header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = f"@{instrument}:{run}:{flowcell}:1:1101:19795:1000 {read}:N:0:AATGAGCG"
    with gzip.open(path, "wt") as fh:
        fh.write(f"{header}\nACGTACGT\n+\nIIIIIIII\n")


def build_two_flowcell_case(root):
    """Two flowcells, identical filenames, one file per run_??? dir.

    Dir order deliberately does not track flowcell -- this is what the real
    work dir looked like, and what defeats positional pairing.
    """
    name = "CITE_CD34_NEU_exp3_libC_HTO_S3_L002_{}_001.fastq.gz"
    layout = [
        ("run_005", "R1", "502", "222KF7NNX"),
        ("run_006", "R2", "501", "222KG5FNX"),
        ("run_007", "R1", "501", "222KG5FNX"),
        ("run_008", "R2", "502", "222KF7NNX"),
    ]
    staged = []
    for run_dir, read, run_no, fc in layout:
        p = root / run_dir / name.format(read)
        write_fastq(p, "VH00206", run_no, fc, read[-1])
        staged.append(p)
    return sorted(staged)


def old_pair(staged):
    """Pre-fix logic: two independent sorts zipped by position."""
    import re
    return list(zip(sorted(p for p in staged if re.search(r"_R1_", p.name)),
                    sorted(p for p in staged if re.search(r"_R2_", p.name))))


def test_pairs_within_flowcell(root):
    staged = build_two_flowcell_case(root)
    pairs, unmatched = pair_reads(staged)
    assert not unmatched, f"unexpected unmatched reads: {unmatched}"
    assert len(pairs) == 2, f"expected 2 pairs, got {len(pairs)}"
    for r1, r2 in pairs:
        assert flowcell_of(r1) == flowcell_of(r2), (
            f"cross-flowcell pair: {r1.parent.name}/{r1.name} {flowcell_of(r1)} "
            f"+ {r2.parent.name}/{r2.name} {flowcell_of(r2)}")


def test_old_logic_actually_mispaired(root):
    """Guards the test itself: the fixture must defeat positional pairing."""
    staged = build_two_flowcell_case(root)
    assert any(flowcell_of(a) != flowcell_of(b) for a, b in old_pair(staged)), (
        "fixture does not reproduce the bug; positional zip paired it correctly")


def test_unmatched_r1_is_reported(root):
    staged = build_two_flowcell_case(root)
    lone = [p for p in staged if "_R1_" in p.name][:1]
    pairs, unmatched = pair_reads(lone)
    assert pairs == [] and len(unmatched) == 1, (
        f"lone R1 should be reported unmatched, got {pairs} {unmatched}")


def test_lane_numbering_and_rename(root):
    staged = build_two_flowcell_case(root)
    final = root / "fastq_all" / "hto"
    lanes, sample_id = stage_modality(staged, final, min_reads=1,
                                      counter=lambda p, m: m)
    assert lanes == 2, f"expected 2 lanes, got {lanes}"
    assert sample_id == "CITE_CD34_NEU_exp3_libC_HTO", sample_id

    written = sorted(p.name for p in final.iterdir())
    assert written == [
        f"{sample_id}_S1_L001_R1_001.fastq.gz",
        f"{sample_id}_S1_L001_R2_001.fastq.gz",
        f"{sample_id}_S1_L002_R1_001.fastq.gz",
        f"{sample_id}_S1_L002_R2_001.fastq.gz",
    ], written

    # The decisive check: each written lane must be internally consistent.
    for lane in ("L001", "L002"):
        r1 = final / f"{sample_id}_S1_{lane}_R1_001.fastq.gz"
        r2 = final / f"{sample_id}_S1_{lane}_R2_001.fastq.gz"
        assert flowcell_of(r1) == flowcell_of(r2), (
            f"{lane} mixes flowcells: {flowcell_of(r1)} vs {flowcell_of(r2)}")


def test_placeholder_pairs_skipped(root):
    staged = build_two_flowcell_case(root)
    final = root / "fastq_all" / "skipme"
    lanes, sample_id = stage_modality(staged, final, min_reads=10000,
                                      counter=lambda p, m: 0)
    assert lanes == 0 and sample_id is None, (lanes, sample_id)
    assert not list(final.iterdir()), "placeholder reads should not be staged"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        root = Path(tempfile.mkdtemp(prefix="oscar_stage_"))
        # stage_modality logs skipped placeholders to stderr; that is expected
        # noise, so buffer it and only reveal it if the test actually fails.
        captured, sys.stderr = sys.stderr, io.StringIO()
        try:
            fn(root)
        except Exception:
            noise, sys.stderr = sys.stderr.getvalue(), captured
            if noise:
                print(noise, file=sys.stderr, end="")
            raise
        else:
            sys.stderr = captured
        finally:
            shutil.rmtree(root, ignore_errors=True)
    print(f"OK: all stage_fastqs self-checks passed ({len(tests)} tests)")


if __name__ == "__main__":
    main()
