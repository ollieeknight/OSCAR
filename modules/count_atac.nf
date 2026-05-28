// ─── CELLRANGER_ATAC ──────────────────────────────────────────────────────────
// Handles: DOGMA-ATAC, Multiome-ATAC, standalone ATAC, ASAP-ATAC.
// Output dir: {library_id}_ATAC  (mirrors current OSCAR naming)
process CELLRANGER_ATAC {
    tag "$meta.library_id"
    label 'process_high'
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

has_r3     = any("_R3_" in p.name for p in staged)
chunk_size = 3 if has_r3 else 2

if len(staged) % chunk_size != 0:
    print(f"[cellranger_atac] ERROR: {len(staged)} files not divisible by chunk_size={chunk_size}", file=sys.stderr)
    sys.exit(1)

final_dir = Path("fastq_all/atac")
final_dir.mkdir(parents=True, exist_ok=True)
lane = 0

for i in range(0, len(staged), chunk_size):
    chunk = staged[i:i + chunk_size]
    r1, r2 = chunk[0], chunk[1]
    r3     = chunk[2] if has_r3 else None
    assert "_R1_" in r1.name, f"Expected R1: {r1}"
    assert "_R2_" in r2.name, f"Expected R2: {r2}"
    if r3:
        assert "_R3_" in r3.name, f"Expected R3: {r3}"
    n = count_reads(r1)
    if n < min_reads:
        print(f"[cellranger_atac] skip {r1.parent.name}: {n} reads", file=sys.stderr)
        continue
    lane += 1
    r1.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R1_001.fastq.gz")
    r2.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R2_001.fastq.gz")
    if r3:
        r3.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R3_001.fastq.gz")

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
