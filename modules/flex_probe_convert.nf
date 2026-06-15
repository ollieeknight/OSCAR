// ─── FLEX_PROBE_PREPARE ────────────────────────────────────────────────────────
// Merges standard 10x probe CSV (comma-delimited) with custom probe CSV
// (semicolon-delimited, OSCAR format) and emits both cellranger and cyto formats.
// Runs once per pipeline invocation (local executor, no SLURM job).
process FLEX_PROBE_PREPARE {
    tag "probe_merge"
    label 'process_low'

    input:
    tuple path(standard_probe_csv), path(custom_probe_csv)

    output:
    path "merged_probes_cellranger.csv", emit: probe_csv_cr
    path "merged_probes_cyto.tsv",       emit: probe_tsv_cyto
    path "versions.yml",                 emit: versions

    script:
    """
    python3 << 'PYEOF'
import csv
from pathlib import Path

standard_rows = []
std_path = Path("${standard_probe_csv}")
if std_path.name != "NO_FILE":
    with open(std_path) as f:
        reader = csv.DictReader(f)
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

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    VERSIONS_EOF
    """
}

// ─── FLEX_BARCODE_EXTRACT ─────────────────────────────────────────────────────
// Extracts the probe barcode reference file from the bundled cellranger container.
// Auto-selects v1 (BC001) or v2 (A-A01) file based on chemistry. Runs once per
// pipeline invocation; result is Nextflow-cached across subsequent runs.
process FLEX_BARCODE_EXTRACT {
    tag "barcode_extract"
    label 'process_low'
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "probe_barcodes.txt", emit: barcodes
    path "versions.yml",       emit: versions

    script:
    def bc_file = chemistry ==~ /Flex-v2.*/
        ? "flex-v2-384.txt"
        : "probe-barcodes-fixed-rna-profiling-rna-r1.txt"
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/translation/${bc_file} probe_barcodes.txt

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oP '(?<=cellranger )[0-9.]+' || echo "unknown")
    VERSIONS_EOF
    """
}

// ─── FLEX_WHITELIST_EXTRACT ───────────────────────────────────────────────────
// Extracts the cell barcode whitelist from the bundled cellranger container.
// Flex v1: 737K-flex-v2.txt.gz (confusingly named, but correct for v1 chemistry).
// Flex v2 (GEM-X): 737K-fixed-rna-profiling.txt.gz.
process FLEX_WHITELIST_EXTRACT {
    tag "whitelist_extract"
    label 'process_low'
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "cb_whitelist.txt.gz", emit: whitelist
    path "versions.yml",        emit: versions

    script:
    def wl_file = chemistry ==~ /Flex-v2.*/
        ? "737K-fixed-rna-profiling.txt.gz"
        : "737K-flex-v2.txt.gz"
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/${wl_file} cb_whitelist.txt.gz

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oP '(?<=cellranger )[0-9.]+' || echo "unknown")
    VERSIONS_EOF
    """
}

// ─── FLEX_SAMPLE_PREPARE ──────────────────────────────────────────────────────
// Expands probe barcode IDs (BC001-BC016) from flex_samples_file to the actual
// 8bp variant sequences that cyto needs for sample demultiplexing.
// Runs once per pipeline invocation (local executor, no SLURM job).
process FLEX_SAMPLE_PREPARE {
    tag "sample_barcodes"
    label 'process_low'

    input:
    path samples_file
    path probe_barcodes_ref

    output:
    path "cyto_probe_barcodes.txt", emit: cyto_barcodes
    path "versions.yml",            emit: versions

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

    cat <<-VERSIONS_EOF > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    VERSIONS_EOF
    """
}
