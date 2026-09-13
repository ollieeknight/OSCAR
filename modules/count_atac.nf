// ─── CELLRANGER_ATAC ─────────────────────────────────────────────────────────
// Handles: DOGMA-ATAC, Multiome-ATAC, standalone ATAC, ASAP-ATAC.
// Output dir: {library_id}_ATAC
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
import os, re, shutil, sys
from pathlib import Path

# Nextflow puts bin/ on PATH, not PYTHONPATH.
sys.path.insert(0, os.path.dirname(shutil.which("stage_fastqs.py")))
from stage_fastqs import count_reads, flowcell_of

sample_id = "${meta.id}"
min_reads  = ${min_reads}

staged = sorted(
    p for p in Path(".").glob("fastqs/atac/run_*/*")
    if p.name != "NO_FILE"
)

# Group the R1/R2/R3 of one lane by flowcell and filename, not by staged dir:
# stageAs gives every file its own run_??? dir, so the dir name records
# input-list position rather than which flowcell a read came from.
runs = {}
for p in staged:
    m = re.search(r'_(R[123])_', p.name)
    if not m:
        continue                       # index reads (I1/I2) are not cellranger-atac inputs
    key = (flowcell_of(p), re.sub(r'_R[123]_', '_R_', p.name))
    runs.setdefault(key, {})[m.group(1)] = p

if not runs:
    print("[cellranger_atac] ERROR: no R1/R2/R3 FASTQs staged", file=sys.stderr)
    sys.exit(1)

final_dir = Path("fastq_all/atac")
final_dir.mkdir(parents=True, exist_ok=True)
lane = 0

for key in sorted(runs):
    reads    = runs[key]
    fq_group = reads["R1"].name
    # cellranger-atac needs all three reads per lane: R1 (insert), R2 (cell
    # barcode) and R3 (insert mate). Checking only R1/R2 let a flowcell that
    # demuxed without R3 through, writing a lane whose read set differs from
    # every other lane's.
    missing  = [t for t in ("R1", "R2", "R3") if t not in reads]
    if missing:
        print(f"[cellranger_atac] ERROR: {fq_group} missing {missing}", file=sys.stderr)
        sys.exit(1)
    n = count_reads(reads["R1"], min_reads)
    if n < min_reads:
        print(f"[cellranger_atac] skip {fq_group}: {n} reads", file=sys.stderr)
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
