#!/usr/bin/env bash
# check_undetermined_barcodes.sh
# Example:
#   bash check_undetermined_barcodes.sh \
#     unknown_bcl \
#     unknown_bcl_indices

set -euo pipefail

BCL_DIR="${1:?Usage: $0 <BCL_DIR> <OUT_DIR>}"
OUT_DIR="${2:?Usage: $0 <BCL_DIR> <OUT_DIR>}"
IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-nf-core-bclconvert-4.4.6.img"

CYCLES=$(python3 - "$BCL_DIR/RunInfo.xml" << 'PYEOF'
import sys, xml.etree.ElementTree as ET
tree = ET.parse(sys.argv[1])
reads = tree.findall('.//Read')
result = []
for r in reads:
    nc = r.get('NumCycles')
    is_index = r.get('IsIndexedRead', 'N')
    result.append((int(nc), is_index == 'Y'))
print(','.join(f"{nc},{int(idx)}" for nc, idx in result))
PYEOF
)

OVERRIDE=$(python3 - "$CYCLES" << 'PYEOF'
import sys
parts = sys.argv[1].split(',')
reads = [(int(parts[i]), parts[i+1]=='1') for i in range(0, len(parts), 2)]
mask = []
for nc, is_index in reads:
    if is_index:
        mask.append(f"I{nc}")
    else:
        mask.append(f"N{nc}")
print(';'.join(mask))
PYEOF
)

mkdir -p "$OUT_DIR"
cat > "$OUT_DIR/SampleSheet.csv" << EOF
[Header]
FileFormatVersion,2

[BCLConvert_Settings]
SoftwareVersion,4.4.6
OverrideCycles,${OVERRIDE}

[BCLConvert_Data]
Sample_ID,index,index2
DUMMY_NOMATCH,AAAAAAAAAA,AAAAAAAAAA
EOF

apptainer exec \
    --bind "$BCL_DIR" \
    --bind "$OUT_DIR" \
    "$IMG" \
    bcl-convert \
        --bcl-input-directory  "$BCL_DIR" \
        --output-directory     "$OUT_DIR/fastq" \
        --sample-sheet         "$OUT_DIR/SampleSheet.csv" \
        --first-tile-only      true \
        --bcl-only-matched-reads false \
        --force

REPORT="$OUT_DIR/fastq/Reports/Top_Unknown_Barcodes.csv"
if [[ -f "$REPORT" ]]; then
    # Show top 20, formatted
    python3 - "$REPORT" << 'PYEOF'
import sys, csv
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f))
rows.sort(key=lambda r: -int(r.get('Count','0').replace(',','')))
print(f"{'index':<12} {'index2':<12} {'Count':>10}  {'Sample'}")
print('-' * 60)
for r in rows[:20]:
    print(f"{r.get('index',''):<12} {r.get('index2',''):<12} {r.get('Count',''):>10}  {r.get('Sample','')}")
PYEOF
else
    echo "ERROR: report not found at $REPORT"
fi
