#!/usr/bin/env bash
# All OSCAR self-checks. Run: bash tests/run_all.sh
# No cluster, no containers, no sequencing data required.
set -euo pipefail
cd "$(dirname "$0")/.."

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "── lib/ helpers ──────────────────────────────────────────"
nextflow run tests/test_lib.nf -w "$WORK" | grep -E "^OK:|ERROR"

echo "── FASTQ staging ─────────────────────────────────────────"
python3 tests/test_stage_fastqs.py

echo "── publish layout ────────────────────────────────────────"
python3 tests/test_publish_layout.py

echo "── pipeline parses ───────────────────────────────────────"
# --help exits non-zero by design (it stops at parameter validation), so
# capture output first rather than letting `set -e` abort on the exit code.
help_out=$(nextflow run main.nf --help -w "$WORK" 2>&1 || true)
if grep -q "samplesheet is required" <<<"$help_out"; then
    echo "OK: main.nf compiles and reaches parameter validation"
else
    echo "FAIL: main.nf did not compile" >&2
    echo "$help_out" | tail -20 >&2
    exit 1
fi

echo
echo "All self-checks passed."
