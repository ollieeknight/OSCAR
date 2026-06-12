// ─── CELLRANGER_MULTI ─────────────────────────────────────────────────────────
// Merges former MULTI_CONFIG + CELLRANGER_MULTI into one process.
// Inputs: 6 path() inputs staged per-file into fastqs/{mod}/run_???/
// Python: read-count filters placeholder FASTQs, renames to L001/L002/...,
//         writes multi_config.csv, then cellranger multi runs.
process CELLRANGER_MULTI {
    tag "$library_id"
    label 'process_high'
    container "${params.container_cellranger}"
    publishDir { "${params.outdir}/${metas[0].run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas), val(config_header), path(adt_csv),
          val(flex_samples_content),
          path(gex_fastqs,    stageAs: "fastqs/gex/run_???/*"),
          path(adt_fastqs,    stageAs: "fastqs/adt/run_???/*"),
          path(hto_fastqs,    stageAs: "fastqs/hto/run_???/*"),
          path(vdj_t_fastqs,  stageAs: "fastqs/vdj_t/run_???/*"),
          path(vdj_b_fastqs,  stageAs: "fastqs/vdj_b/run_???/*"),
          path(crispr_fastqs, stageAs: "fastqs/crispr/run_???/*")

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs

    script:
    def min_reads = 10000
    """
    cat > config_header.txt << 'OSCAR_CR_HEADER_EOF'
${config_header}
OSCAR_CR_HEADER_EOF

    cat > flex_samples.txt << 'OSCAR_FLEX_SAMPLES_EOF'
${flex_samples_content}
OSCAR_FLEX_SAMPLES_EOF

    python3 << 'OSCAR_PYEOF'
import re, subprocess, sys
from pathlib import Path

library_id = "${library_id}"
min_reads  = ${min_reads}

with open("config_header.txt") as fh:
    config_header = fh.read().strip()

with open("flex_samples.txt") as fh:
    flex_samples_content = fh.read().strip()

MODALITY_FEATURE = {
    "gex":    "Gene Expression",
    "adt":    "Antibody Capture",
    "hto":    "Antibody Capture",
    "vdj_t":  "VDJ-T",
    "vdj_b":  "VDJ-B",
    "crispr": "CRISPR Guide Capture",
}

def count_reads(r1):
    result = subprocess.run(
        f"zcat {r1} | head -n {min_reads * 4} | awk 'NR%4==1' | wc -l",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip() or "0")

work_dir  = str(Path.cwd())
lib_lines = []

for mod, feat in MODALITY_FEATURE.items():
    staged = sorted(
        p for p in Path(".").glob(f"fastqs/{mod}/run_*/*")
        if p.name != "NO_FILE"
    )
    if not staged:
        continue

    r1_files = sorted(p for p in staged if re.search(r'_R1_', p.name))
    r2_files = sorted(p for p in staged if re.search(r'_R2_', p.name))
    if len(r1_files) != len(r2_files):
        print(f"[cellranger_multi] ERROR: R1/R2 count mismatch for {mod}: {len(r1_files)} R1s vs {len(r2_files)} R2s", file=sys.stderr)
        sys.exit(1)
    if not r1_files:
        continue

    final_dir = Path("fastq_all") / mod
    final_dir.mkdir(parents=True, exist_ok=True)
    lane      = 0
    sample_id = None

    for r1, r2 in zip(r1_files, r2_files):
        n = count_reads(r1)
        if n < min_reads:
            print(f"[cellranger_multi] skip {mod}/{r1.parent.name}: {n} reads", file=sys.stderr)
            continue
        lane += 1
        sample_id = re.sub(r'_S\\d+.*', '', r1.name)
        r1.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R1_001.fastq.gz")
        r2.rename(final_dir / f"{sample_id}_S1_L{lane:03d}_R2_001.fastq.gz")

    if lane > 0 and sample_id:
        lib_lines.append(f"{sample_id},{work_dir}/fastq_all/{mod},{feat}")

cfg = config_header
if lib_lines:
    cfg += "\\n\\n[libraries]\\nfastq_id,fastqs,feature_types\\n"
    cfg += "\\n".join(lib_lines)
if flex_samples_content:
    cfg += "\\n" + flex_samples_content
cfg = re.sub(r'\\n{3,}', '\\n\\n', cfg)

with open("multi_config.csv", "w") as fh:
    fh.write(cfg + "\\n")

print("[cellranger_multi] Config written:", file=sys.stderr)
with open("multi_config.csv") as fh:
    print(fh.read(), file=sys.stderr)
OSCAR_PYEOF

    cellranger multi \\
        --id        "${library_id}" \\
        --csv       multi_config.csv \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}

    rm -rf "${library_id}/SC_MULTI_CS" "${library_id}/_"*
    """
}
