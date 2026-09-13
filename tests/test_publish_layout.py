#!/usr/bin/env python3
"""Checks the publish-layout invariant. Run: python3 tests/test_publish_layout.py

Per sequencing folder there are only ever two generated directories:

    {run}_fastq   always -- FASTQs plus all read-level QC
    {run}_outs    only for the final run of a resequencing chain

Every published file must land inside one of those. Regression target: FALCO
and MULTIQC published to ${params.outdir}/${run_name}_fastq, so merging four
flowcells scattered R462_fastq/, R463_fastq/, R501_fastq/ into the R502 folder.
"""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MODULES = sorted((REPO / "modules").glob("*.nf"))

# publishDir { "..." } / publishDir "..." -- capture the path expression.
PUBLISH = re.compile(r"publishDir\s*(?:\{\s*)?\"([^\"]+)\"")


def publish_targets():
    """[(module_name, line_no, path_expression)] for every publishDir."""
    found = []
    for mod in MODULES:
        for i, line in enumerate(mod.read_text().splitlines(), 1):
            m = PUBLISH.search(line)
            if m:
                found.append((mod.name, i, m.group(1)))
    return found


def test_every_publish_is_fastq_or_outs():
    """Each target's first path segment must be a _fastq or _outs directory."""
    bad = []
    for mod, line, expr in publish_targets():
        # Root is whatever precedes the first generated directory component.
        if "_outs" in expr or "_fastq" in expr:
            continue
        # A target derived from an already-validated dir is fine.
        if "${fastq_dir}" in expr or "${bcl_parent}" in expr:
            continue
        bad.append(f"{mod}:{line} -> {expr}")
    assert not bad, (
        "publishDir targets outside {run}_fastq / {run}_outs:\n  " + "\n  ".join(bad))


def process_body(name):
    """Body of a process up to its script block, found in whichever module holds it."""
    for mod in MODULES:
        block = mod.read_text().split(f"process {name} {{", 1)
        if len(block) == 2:
            return block[1].split("\n    script:", 1)[0]
    raise AssertionError(f"process {name} not found in modules/")


def test_fastq_qc_publishes_beside_its_own_fastqs():
    """Read-level QC must derive its path from fastq_dir, not params.outdir.

    Deriving from params.outdir is what created the stray {run}_fastq trees:
    with --extra_bcl_dirs every run's QC landed under the primary run's folder.
    """
    offenders = []
    for proc in ("FASTP", "MULTIQC"):
        body = process_body(proc)
        m = PUBLISH.search(body)
        assert m, f"{proc} has no publishDir"
        if "${fastq_dir}" not in m.group(1):
            offenders.append(f"{proc} -> {m.group(1)}")
    assert not offenders, (
        "read-level QC must publish under ${fastq_dir}:\n  " + "\n  ".join(offenders))


def test_bclconvert_publishes_beside_source_flowcell():
    body = process_body("BCLCONVERT")
    assert "${bcl_parent}/${run}_fastq" in body, (
        "BCLCONVERT must publish to {bcl_parent}/{run}_fastq so each flowcell's "
        "reads land beside their own BCL directory")


def test_no_falco_references_remain():
    """Falco was replaced by fastp in report-only mode."""
    stale = []
    for f in list(REPO.glob("*.nf")) + list(REPO.glob("*.config")) + \
             list((REPO / "modules").glob("*.nf")) + \
             list((REPO / "subworkflows").glob("*.nf")) + \
             list((REPO / "assets").glob("*.yaml")):
        if re.search(r"falco", f.read_text(), re.I):
            stale.append(f.relative_to(REPO).as_posix())
    assert not stale, f"stale falco references in: {stale}"


def test_outdir_is_derived_not_hardcoded():
    """--outdir must default to the sibling of the run's input folder.

    Rooting on a fixed default put R462_fastq/, R463_fastq/ and R501_fastq/
    inside 2608_R502_Sergio when four flowcells were merged.
    """
    cfg = (REPO / "nextflow.config").read_text()
    block = re.search(r"outdir\s*=\s*(.+?)\n\s*run_from", cfg, re.S)
    assert block, "could not locate the outdir parameter in nextflow.config"
    expr = block.group(1)
    for anchor in ("params.bcl_dir", "params.fastq_dir", "params.outs_dir"):
        assert anchor in expr, (
            f"outdir must derive from {anchor} so each entry point roots "
            f"output beside its own input folder; got: {expr.strip()}")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"OK: all publish-layout self-checks passed ({len(tests)} tests)")


if __name__ == "__main__":
    main()
