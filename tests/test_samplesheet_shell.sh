#!/usr/bin/env bash
# Self-checks for the shell half of GENERATE_SAMPLESHEET: RunInfo.xml cycle
# parsing, OverrideCycles wildcard expansion, and the global-vs-per-sample
# layout of the emitted SampleSheet.csv.
#
# The script under test is extracted from modules/demux.nf rather than copied,
# so this cannot silently drift from the process it is meant to cover.
set -euo pipefail
cd "$(dirname "$0")/.."

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ── Extract the process script body and substitute the Nextflow interpolations
# the Groovy half would have filled in. This mirrors a mixed CITE group:
# two 10bp DI GEX rows and two 8bp SI ADT/HTO rows.
python3 - "$TMP" <<'PYEOF'
import re, sys, os
tmp = sys.argv[1]
src = open('modules/demux.nf').read()
block = src[src.index('// ─── GENERATE_SAMPLESHEET ─'):src.index('// ─── CLEAN_FASTQ_DIR ─')]
# The script body is the last triple-quoted string in the process.
body = block[block.index('"""') + 3:block.rindex('"""')]

specs = "\n".join([
    "CITE_A_GEX\tGTGGATCGCA\tACACTAAGGC\tY28N*;I10N*;I10N*;Y90N*\tY28N*;I10N*;Y90N*",
    "CITE_B_GEX\tACAGACCTTT\tGTATTCCACC\tY28N*;I10N*;I10N*;Y90N*\tY28N*;I10N*;Y90N*",
    "CITE_A_ADT\tATTCAGAA\t-\tY28N*;I8N2;N*;Y90N*\tY28N*;I8N*;Y90N*",
    "CITE_A_HTO\tTTCTGAAT\t-\tY28N*;I8N2;N*;Y90N*\tY28N*;I8N*;Y90N*",
])

# Undo Nextflow's shell escaping the way Nextflow itself does: a single
# left-to-right pass over the escape pairs. Doing it as two sequential
# str.replace calls corrupts \\t (the backslash is consumed by the first pass).
def unescape(text):
    out, i = [], 0
    while i < len(text):
        if text[i] == '\\' and i + 1 < len(text):
            out.append(text[i + 1]); i += 2
        else:
            out.append(text[i]); i += 1
    return ''.join(out)

body = unescape(body)
body = body.replace('${bcl_dir}', os.path.join(tmp, 'bcl'))
body = body.replace('${spec_tsv}', specs)
body = body.replace('${per_sample}', 'true')
body = body.replace('${is_dual}', 'true')
open(os.path.join(tmp, 'mixed.sh'), 'w').write(body)

# Uniform group: all rows share one mask, so no per-sample column.
uni = "\n".join([
    "CITE_A_GEX\tGTGGATCGCA\tACACTAAGGC\tY28N*;I10N*;I10N*;Y90N*\tY28N*;I10N*;Y90N*",
    "CITE_B_GEX\tACAGACCTTT\tGTATTCCACC\tY28N*;I10N*;I10N*;Y90N*\tY28N*;I10N*;Y90N*",
])
body2 = block[block.index('"""') + 3:block.rindex('"""')]
body2 = unescape(body2)
body2 = body2.replace('${bcl_dir}', os.path.join(tmp, 'bcl'))
body2 = body2.replace('${spec_tsv}', uni)
body2 = body2.replace('${per_sample}', 'false')
body2 = body2.replace('${is_dual}', 'true')
open(os.path.join(tmp, 'uniform.sh'), 'w').write(body2)
PYEOF

# ── A 4-read NovaSeq X style RunInfo.xml: 28 / 10 / 10 / 90.
mkdir -p "$TMP/bcl"
cat > "$TMP/bcl/RunInfo.xml" <<'EOF'
<?xml version="1.0"?>
<RunInfo><Run><Reads>
  <Read Number="1" NumCycles="28" IsIndexedRead="N"/>
  <Read Number="2" NumCycles="10" IsIndexedRead="Y"/>
  <Read Number="3" NumCycles="10" IsIndexedRead="Y"/>
  <Read Number="4" NumCycles="90" IsIndexedRead="N"/>
</Reads></Run></RunInfo>
EOF

run_case() {
    local name="$1"
    local dir="$TMP/run_$name"
    mkdir -p "$dir"
    ( cd "$dir" && bash "$TMP/$name.sh" >/dev/null )
    cat "$dir/SampleSheet.csv"
}

# ── mixed group ───────────────────────────────────────────────────────────
MIXED=$(run_case mixed)

grep -q '^Read1Cycles,28$'  <<<"$MIXED" || { echo "FAIL: Read1Cycles"; exit 1; }
grep -q '^Index1Cycles,10$' <<<"$MIXED" || { echo "FAIL: Index1Cycles"; exit 1; }
grep -q '^Index2Cycles,10$' <<<"$MIXED" || { echo "FAIL: Index2Cycles"; exit 1; }
grep -q '^Read2Cycles,90$'  <<<"$MIXED" || { echo "FAIL: Read2Cycles"; exit 1; }
echo "OK: cycle counts read from RunInfo.xml"

# Illumina: a setting is global OR per-sample, never both.
if grep -q '^OverrideCycles,' <<<"$MIXED"; then
    echo "FAIL: mixed sheet has a global OverrideCycles alongside the per-sample column"
    exit 1
fi
grep -q '^Sample_ID,Index,Index2,OverrideCycles$' <<<"$MIXED" || {
    echo "FAIL: expected per-sample OverrideCycles column header"; echo "$MIXED"; exit 1; }
# index2 is not used for demultiplexing on the SI rows, so a global
# BarcodeMismatchesIndex2 would not be valid for this sheet.
if grep -q '^BarcodeMismatchesIndex2,' <<<"$MIXED"; then
    echo "FAIL: mixed sheet sets BarcodeMismatchesIndex2 while SI rows mask i5"
    exit 1
fi
grep -q '^BarcodeMismatchesIndex1,1$' <<<"$MIXED" || {
    echo "FAIL: mixed sheet lost BarcodeMismatchesIndex1"; exit 1; }
echo "OK: mixed group emits per-sample column and no global setting"

# Wildcards must be gone — bcl-convert 4.x rejects '*'.
if grep -q '\*' <<<"$MIXED"; then
    echo "FAIL: unexpanded wildcard left in sheet"; echo "$MIXED"; exit 1
fi
echo "OK: all OverrideCycles wildcards expanded"

# The whole point: each row's declared index cycles match its own index length.
python3 - <<PYEOF
import re, sys
sheet = """$MIXED"""
data = sheet.split('[BCLConvert_Data]')[1].strip().splitlines()
hdr, rows = data[0].split(','), [r.split(',') for r in data[1:] if r.strip()]
assert hdr == ['Sample_ID', 'Index', 'Index2', 'OverrideCycles'], hdr
seen = {}
for r in rows:
    sid, i7, i5, mask = r[0], r[1], r[2], r[3]
    parts = mask.split(';')
    i1 = int(re.match(r'I(\d+)', parts[1]).group(1))
    assert i1 == len(i7), f"{sid}: mask {mask} declares I{i1} but index {i7} is {len(i7)}bp"
    # Read1/Read2 keep their full length on every row.
    assert parts[0] == 'Y28', f"{sid}: read1 mask {parts[0]}"
    assert parts[3] == 'Y90', f"{sid}: read2 mask {parts[3]}"
    if i5:
        # Dual-index row: i5 cycles declared must match the i5 length.
        i2 = int(re.match(r'I(\d+)', parts[2]).group(1))
        assert i2 == len(i5), f"{sid}: mask {mask} vs index2 {i5}"
    else:
        # Single-index row: i5 fully masked, never demultiplexed on.
        assert parts[2] == 'N10', f"{sid}: expected i5 fully masked, got {parts[2]}"
    seen[sid] = mask

assert seen['CITE_A_GEX'] == 'Y28;I10;I10;Y90', seen['CITE_A_GEX']
assert seen['CITE_A_ADT'] == 'Y28;I8N2;N10;Y90', seen['CITE_A_ADT']
assert seen['CITE_A_HTO'] == 'Y28;I8N2;N10;Y90', seen['CITE_A_HTO']

# Every row of a lane must consume the same total cycles per read.
totals = set()
for r in rows:
    tot = tuple(sum(int(n) for n in re.findall(r'\d+', p)) for p in r[3].split(';'))
    totals.add(tot)
assert len(totals) == 1, f"rows disagree on total cycles per read: {totals}"
assert totals.pop() == (28, 10, 10, 90)
print("OK: per-row masks valid, i5 masked on SI rows, all rows consume 28/10/10/90")
PYEOF

# ── uniform group keeps the global setting and no column ──────────────────
UNIFORM=$(run_case uniform)
grep -q '^OverrideCycles,Y28;I10;I10;Y90$' <<<"$UNIFORM" || {
    echo "FAIL: uniform sheet lost its global OverrideCycles"; echo "$UNIFORM"; exit 1; }
grep -q '^Sample_ID,Index,Index2$' <<<"$UNIFORM" || {
    echo "FAIL: uniform sheet should not have a per-sample column"; echo "$UNIFORM"; exit 1; }
grep -q '^BarcodeMismatchesIndex2,1$' <<<"$UNIFORM" || {
    echo "FAIL: uniform dual-index sheet should keep BarcodeMismatchesIndex2"; exit 1; }
echo "OK: uniform group keeps global OverrideCycles, no per-sample column"

echo "OK: all samplesheet shell self-checks passed"
