// ─── CYTO_FLEX ────────────────────────────────────────────────────────────────
// Maps Flex GEX reads to probe sequences using cyto (ARC Institute).
// Produces per-sample MTX outputs (one dir per probe barcode) for probe-level QC.
// Does NOT require a reference genome — maps reads directly to probe sequences.
process CYTO_FLEX {
    tag "$library_id"
    label 'process_high'
    container "${params.container_cyto}"
    publishDir { "${params.outdir}/${metas[0].run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas),
          path(probe_tsv_cyto),
          path(cyto_probe_barcodes),
          path(cb_whitelist),
          path(gex_fastqs, stageAs: "fastqs/gex/run_???/*")

    output:
    tuple val(library_id), val(metas), path("${library_id}_cyto"), emit: counts
    path "versions.yml", emit: versions

    script:
    def min_reads = 10000
    """
    python3 << 'OSCAR_PYEOF'
import re, subprocess, sys
from pathlib import Path

min_reads = ${min_reads}
staged    = sorted(p for p in Path(".").glob("fastqs/gex/run_*/*") if p.name != "NO_FILE")
r1_files  = sorted(p for p in staged if re.search(r'_R1_', p.name))
r2_files  = sorted(p for p in staged if re.search(r'_R2_', p.name))

def count_reads(r1):
    result = subprocess.run(
        f"zcat {r1} | head -n {min_reads * 4} | awk 'NR%4==1' | wc -l",
        shell=True, capture_output=True, text=True
    )
    return int(result.stdout.strip() or "0")

pairs = []
for r1, r2 in zip(r1_files, r2_files):
    n = count_reads(r1)
    if n >= min_reads:
        pairs.extend([str(r1), str(r2)])
    else:
        print(f"[cyto_flex] skip {r1.parent.name}: {n} reads (<{min_reads})", file=sys.stderr)

if not pairs:
    print("[cyto_flex] ERROR: no FASTQ pairs passed read threshold", file=sys.stderr)
    sys.exit(1)

with open("fastq_pairs.txt", "w") as f:
    f.write(" ".join(pairs))
OSCAR_PYEOF

    FASTQ_PAIRS=\$(cat fastq_pairs.txt)

    PROBES_ARG=""
    if [ "${cyto_probe_barcodes}" != "NO_FILE" ]; then
        PROBES_ARG="-p ${cyto_probe_barcodes}"
    fi

    cyto workflow gex \\
        -c ${probe_tsv_cyto} \\
        \${PROBES_ARG} \\
        -w ${cb_whitelist} \\
        --preset gex-v2 \\
        -o ${library_id}_cyto \\
        -F mtx \\
        --no-filter \\
        --memory-limit ${params.flex_cyto_memory_limit} \\
        -T ${task.cpus} \\
        -f \\
        \${FASTQ_PAIRS}

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        cyto: \$(cyto --version 2>&1 | awk '{print \$NF}')
    VERSIONS_EOF
    """
}

// ─── CYTO_RENAME_SAMPLES ──────────────────────────────────────────────────────
// Renames cyto output directories from probe barcode IDs (BC001, BC002, ...)
// to sample names from flex_samples_file (P15_PBMC, H13_PBMC, ...).
// No-op when samples_file is NO_FILE (singleplex / uniplexed run).
process CYTO_RENAME_SAMPLES {
    tag "$library_id"
    label 'process_low'

    input:
    tuple val(library_id), val(metas), path(cyto_out), path(samples_file)

    output:
    tuple val(library_id), val(metas), path(cyto_out), emit: counts
    path "versions.yml", emit: versions

    script:
    """
    python3 << 'PYEOF'
import csv, shutil
from pathlib import Path

samples_path = Path("${samples_file}")
cyto_out     = Path("${cyto_out}")
counts_dir   = cyto_out / "counts"

if samples_path.name == "NO_FILE" or not counts_dir.exists():
    print("[cyto_rename] no samples file or no counts dir — skipping rename")
else:
    # Build per-BC mapping: each individual BC ID → sample_id.
    # probe_barcode_ids may be multi-valued (e.g. "BC001|BC002") for pooled samples;
    # cyto creates one output dir per individual BC, so we expand those here.
    mapping = {}
    with open(samples_path) as f:
        for row in csv.DictReader(f):
            sample = row['sample_id'].strip()
            for bc in row['probe_barcode_ids'].strip().replace('|', ',').split(','):
                mapping[bc.strip()] = sample

    for bc_dir in sorted(counts_dir.iterdir()):
        if bc_dir.name in mapping:
            target = bc_dir.parent / mapping[bc_dir.name]
            shutil.move(str(bc_dir), str(target))
            print(f"[cyto_rename] {bc_dir.name} → {mapping[bc_dir.name]}")
        else:
            print(f"[cyto_rename] no mapping for {bc_dir.name} — left as-is")
PYEOF

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    VERSIONS_EOF
    """
}
