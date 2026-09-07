// ─── CELLRANGER_ATAC ──────────────────────────────────────────────────────────
// Handles: DOGMA-ATAC, Multiome-ATAC, standalone ATAC, ASAP-ATAC.
// Output dir: {library_id}_ATAC  (mirrors current OSCAR naming)
process CELLRANGER_ATAC {
    tag "$meta.library_id"
    container "${params.container_cellranger_atac}"
    publishDir { "${params.outdir}/${meta.run_name}_outs" }, mode: 'copy'

    input:
    tuple val(meta), path(atac_fastqs, stageAs: "fastqs/atac/run_???/*")

    output:
    tuple val(meta), path("${meta.library_id}_ATAC/outs"), emit: outs

    script:
    def min_reads  = 10000
    def reference  = meta.species == 'human' ? params.ref_human : params.ref_mouse
    def extra_args = (meta.assay == 'DOGMA') ? "\\\n        --chemistry ARC-v1" : ''
    """
    python3 << 'OSCAR_ATAC_PYEOF'
import re, subprocess, sys
from pathlib import Path

sample_id = "${meta.id}"
min_reads  = ${min_reads}

def count_reads(r1):
    result = subprocess.run(
        f"zcat {r1} | head -n {min_reads * 4} | awk 'NR%4==1' | wc -l",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip() or "0")

staged = sorted(
    p for p in Path(".").glob("fastqs/atac/run_*/*")
    if p.name != "NO_FILE"
)

# Group per staged run_??? dir and key by read tag. Blind positional chunking
# breaks whenever BCL Convert also emits index reads (_I1_/_I2_ sort before _R1_).
runs = {}
for p in staged:
    m = re.search(r'_(R[123])_', p.name)
    if not m:
        continue                       # index reads (I1/I2) are not cellranger-atac inputs
    runs.setdefault(p.parent.name, {})[m.group(1)] = p

if not runs:
    print("[cellranger_atac] ERROR: no R1/R2/R3 FASTQs staged", file=sys.stderr)
    sys.exit(1)

final_dir = Path("fastq_all/atac")
final_dir.mkdir(parents=True, exist_ok=True)
lane = 0

for run_name in sorted(runs):
    reads = runs[run_name]
    missing = [t for t in ("R1", "R2") if t not in reads]
    if missing:
        print(f"[cellranger_atac] ERROR: {run_name} missing {missing}", file=sys.stderr)
        sys.exit(1)
    n = count_reads(reads["R1"])
    if n < min_reads:
        print(f"[cellranger_atac] skip {run_name}: {n} reads", file=sys.stderr)
        continue
    lane += 1
    for tag, p in sorted(reads.items()):
        p.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_{tag}_001.fastq.gz")

if lane == 0:
    print("[cellranger_atac] ERROR: no valid FASTQ lanes found", file=sys.stderr)
    sys.exit(1)

with open("fastq_dir.txt", "w") as fh:
    fh.write(str(Path("fastq_all/atac").resolve()))
OSCAR_ATAC_PYEOF

    cellranger-atac count \\
        --id        "${meta.library_id}_ATAC" \\
        --reference "${reference}" \\
        --fastqs    \$(cat fastq_dir.txt) \\
        --sample    "${meta.id}" \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}${extra_args}

    rm -rf "${meta.library_id}_ATAC/SC_ATAC_COUNTER_CS" "${meta.library_id}_ATAC/_"*
    """
}
