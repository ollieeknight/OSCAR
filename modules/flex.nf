include { get_flex_barcode_file; get_flex_whitelist_file } from '../lib/chemistry'

// ─── Flex (Fixed RNA Profiling) ──────────────────────────────────────────────
// Probe-set preparation and barcode extraction for cellranger multi, plus the
// cyto backend for probe-level QC. Backend chosen by params.flex_backend.

// ─── FLEX_PROBE_PREPARE ──────────────────────────────────────────────────────
// Merges standard 10x probe CSV (comma-delimited) with custom probe CSV
// (semicolon-delimited, OSCAR format) and emits both cellranger and cyto formats.
// Runs once per pipeline invocation (local executor, no SLURM job).
process FLEX_PROBE_PREPARE {
    container "${params.container_python}"

    input:
    tuple path(standard_probe_csv), path(custom_probe_csv)

    output:
    path "merged_probes_cellranger.csv", emit: probe_csv_cr
    path "merged_probes_cyto.tsv",       emit: probe_tsv_cyto

    script:
    """
    python3 << 'PYEOF'
import csv
from pathlib import Path

standard_rows = []
std_path = Path("${standard_probe_csv}")
if std_path.name != "NO_FILE":
    with open(std_path) as f:
        lines = [l for l in f if not l.startswith('#')]
        reader = csv.DictReader(lines)
        for row in reader:
            standard_rows.append(row)

custom_rows = []
cust_path = Path("${custom_probe_csv}")
if cust_path.name != "NO_FILE":
    with open(cust_path) as f:
        reader = csv.DictReader(f, delimiter=';')
        for row in reader:
            custom_rows.append({
                'gene_id':   row['gene_id'],
                'probe_seq': row['probe_seq'],
                'probe_id':  row['probe_id'],
                'included':  'TRUE' if row['included'].upper() == 'TRUE' else 'FALSE',
                'region':    row.get('region', 'spliced'),
                'gene_name': row.get('gene_name', row['gene_id']),
            })

# Custom probes override standard by probe_id
seen = {}
for r in standard_rows + custom_rows:
    seen[r['probe_id']] = r
merged = list(seen.values())

# cellranger CSV: comma-delimited, required cols + region if present
cr_fields = ['gene_id', 'probe_seq', 'probe_id', 'included']
if merged and 'region' in merged[0]:
    cr_fields.append('region')
with open("merged_probes_cellranger.csv", "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=cr_fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(merged)

# cyto GEX library TSV: tab-delimited, NO header, included probes only
# Column order: probe_id\\tgene_name\\tprobe_seq
with open("merged_probes_cyto.tsv", "w") as f:
    for r in merged:
        if r.get('included', 'TRUE').upper() == 'TRUE':
            gene_name = r.get('gene_name', r['gene_id'])
            f.write(f"{r['probe_id']}\\t{gene_name}\\t{r['probe_seq']}\\n")

print(f"Merged {len(merged)} probes ({len(standard_rows)} standard + {len(custom_rows)} custom)")
PYEOF
    """
}

// ─── FLEX_BARCODE_EXTRACT ────────────────────────────────────────────────────
// Extracts the probe barcode reference file from the bundled cellranger container.
// Auto-selects v1 (BC001) or v2 (A-A01) file based on chemistry. Runs once per
// pipeline invocation; result is Nextflow-cached across subsequent runs.
process FLEX_BARCODE_EXTRACT {
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "probe_barcodes.txt", emit: barcodes

    script:
    def bc_file = get_flex_barcode_file(chemistry)
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/translation/${bc_file} probe_barcodes.txt
    """
}

// ─── FLEX_WHITELIST_EXTRACT ──────────────────────────────────────────────────
// Extracts the cell barcode whitelist from the bundled cellranger container.
// Flex v1: 737K-flex-v2.txt.gz (confusingly named, but correct for v1 chemistry).
// Flex v2 (GEM-X): 737K-fixed-rna-profiling.txt.gz.
process FLEX_WHITELIST_EXTRACT {
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "cb_whitelist.txt.gz", emit: whitelist

    script:
    def wl_file = get_flex_whitelist_file(chemistry)
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/${wl_file} cb_whitelist.txt.gz
    """
}

// ─── FLEX_SAMPLE_PREPARE ─────────────────────────────────────────────────────
// Expands probe barcode IDs (BC001-BC016) from flex_samples_file to the actual
// 8bp variant sequences that cyto needs for sample demultiplexing.
// Runs once per pipeline invocation (local executor, no SLURM job).
process FLEX_SAMPLE_PREPARE {
    container "${params.container_python}"

    input:
    path samples_file
    path probe_barcodes_ref

    output:
    path "cyto_probe_barcodes.txt", emit: cyto_barcodes

    script:
    """
    python3 << 'PYEOF'
import csv

# Build BC_ID → [(variant8, canonical8)] from 10x probe barcode reference
# File format: space-separated: variant8 canonical8 barcode_id (no header)
bc_map = {}
with open("${probe_barcodes_ref}") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        variant, canonical, bc_id = parts[0], parts[1], parts[2]
        bc_map.setdefault(bc_id, []).append((variant, canonical))

# Collect used BC IDs from samples file (supports | or , separated multi-BC)
used_bcs = set()
with open("${samples_file}") as f:
    for row in csv.DictReader(f):
        for bc in row['probe_barcode_ids'].strip().replace('|', ',').split(','):
            used_bcs.add(bc.strip())

with open("cyto_probe_barcodes.txt", "w") as out:
    for bc_id in sorted(used_bcs):
        entries = bc_map.get(bc_id, [])
        if not entries:
            print(f"WARNING: no sequences found for {bc_id} in probe barcode reference")
        for variant, canonical in entries:
            out.write(f"{variant}\\t{canonical}\\t{bc_id}\\n")

print(f"Wrote barcode sequences for {len(used_bcs)} probe barcodes")
PYEOF
    """
}

// ─── CYTO_FLEX ───────────────────────────────────────────────────────────────
// Maps Flex GEX reads to probe sequences using cyto (ARC Institute).
// Produces per-sample MTX outputs (one dir per probe barcode) for probe-level QC.
// Does NOT require a reference genome — maps reads directly to probe sequences.
process CYTO_FLEX {
    tag "$library_id"
    container "${params.container_cyto}"
    publishDir { "${params.outdir}/${metas[0].run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas),
          path(probe_tsv_cyto),
          path(cyto_probe_barcodes),
          path(cb_whitelist),
          val(cyto_preset),
          path(gex_fastqs, stageAs: "fastqs/gex/run_???/*")

    output:
    tuple val(library_id), val(metas), path("${library_id}_cyto"), emit: counts

    script:
    def min_reads = 10000
    """
    python3 << 'OSCAR_PYEOF'
import os, shutil, sys
from pathlib import Path

# Nextflow puts bin/ on PATH, not PYTHONPATH.
sys.path.insert(0, os.path.dirname(shutil.which("stage_fastqs.py")))
from stage_fastqs import count_reads, pair_reads

min_reads = ${min_reads}
staged    = sorted(p for p in Path(".").glob("fastqs/gex/run_*/*") if p.name != "NO_FILE")

# Pair on the read header, not on sorted position: two flowcells can stage
# identically-named FASTQs, and zip() would then pair R1 with a foreign R2.
read_pairs, unmatched = pair_reads(staged)
if unmatched:
    print(f"[cyto_flex] ERROR: unmatched reads: {unmatched}", file=sys.stderr)
    sys.exit(1)

pairs = []
for r1, r2 in read_pairs:
    n = count_reads(r1, min_reads)
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
        --preset ${cyto_preset} \\
        -o ${library_id}_cyto \\
        -F mtx \\
        --no-filter \\
        --memory-limit ${params.flex_cyto_memory_limit} \\
        -T ${task.cpus} \\
        -f \\
        \${FASTQ_PAIRS}
    """
}

// ─── CYTO_RENAME_SAMPLES ─────────────────────────────────────────────────────
// Renames cyto output directories from probe barcode IDs (BC001, BC002, ...)
// to sample names from flex_samples_file (P15_PBMC, H13_PBMC, ...).
// No-op when samples_file is NO_FILE (singleplex / uniplexed run).
process CYTO_RENAME_SAMPLES {
    tag "$library_id"

    input:
    tuple val(library_id), val(metas), path(cyto_out), path(samples_file)

    output:
    tuple val(library_id), val(metas), path(cyto_out), emit: counts

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
    """
}
