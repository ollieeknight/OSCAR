// ─── CELLBENDER ───────────────────────────────────────────────────────────────
// Ambient RNA removal. GPU process.
// Source: 05_quality_control.sh:91-139

process CELLBENDER {
    tag "$meta.library_id"
    label 'process_gpu'
    container "${params.container_cellbender}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/cellbender" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("cellbender_out/output.h5"),               emit: h5
    tuple val(meta), path("cellbender_out/output_cell_barcodes.csv"), emit: barcodes
    path "cellbender_out/*"

    script:
    """
    feature_matrix=\$(find -L ${outs_dir} -name 'raw_feature_bc_matrix.h5' | head -1)
    if [ -z "\$feature_matrix" ]; then
        echo "ERROR: raw_feature_bc_matrix.h5 not found in ${outs_dir}" >&2
        exit 1
    fi

    mkdir -p cellbender_out

    cellbender remove-background \\
        --cuda \\
        --input  "\$feature_matrix" \\
        --output cellbender_out/output.h5 \\
        --cpu-threads ${task.cpus} \\
        --checkpoint-mins 10000
    """
}

// ─── CELLSNP_LITE ─────────────────────────────────────────────────────────────
// SNP pileup for donor demultiplexing.
// Two modes: GEX (uses UMI-tagged BAM) and ATAC (no UMI tag).
// Source: 05_quality_control.sh:196-204, 265-273, 350-359

process CELLSNP_LITE {
    tag { "${meta.library_id} (${mode})" }
    label 'process_medium'   // overridden to 32c/64GB/96h via withName
    container "${params.container_cellsnp}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${mode == 'atac' ? "${meta.library_id}_ATAC" : meta.library_id}/vireo" }, mode: 'copy',
               saveAs: { fn -> fn.contains('/') ? file(fn).name : null }

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)
    val(mode)       // 'gex' or 'atac'

    output:
    tuple val(meta), path("cellsnp_${meta.library_id}/"), emit: vcf
    path "cellsnp_${meta.library_id}/*"

    script:
    def out_dir    = (mode == 'atac') ? "${meta.library_id}_ATAC" : meta.library_id
    def umi_flag   = (mode == 'atac') ? '--UMItag None' : ''
    """
    mkdir -p cellsnp_${meta.library_id}

    cellsnp-lite \\
        -s  ${bam} \\
        -b  ${barcodes} \\
        -O  cellsnp_${meta.library_id} \\
        -R  ${params.snp_vcf} \\
        --minMAF   0.1 \\
        --minCOUNT 20 \\
        --gzip \\
        -p  ${task.cpus} \\
        ${umi_flag}
    """
}

// ─── VIREO ────────────────────────────────────────────────────────────────────
// Probabilistic donor demultiplexing.
// Source: 05_quality_control.sh:210-214, 279-283, 365-369

process VIREO {
    tag "$meta.library_id"
    label 'process_medium'   // overridden to 32c/64GB/96h via withName: 'CELLSNP_LITE|VIREO'
    container "${params.container_vireo}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${mode == 'atac' ? "${meta.library_id}_ATAC" : meta.library_id}/vireo" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(cellsnp_dir)
    val(mode)   // 'gex' or 'atac'

    output:
    tuple val(meta), path("vireo_out/donor_ids.tsv"), emit: donor_ids
    path "vireo_out/*"

    script:
    def out_dir = (mode == 'atac') ? "${meta.library_id}_ATAC" : meta.library_id
    """
    mkdir -p vireo_out
    vireo \\
        -c ${cellsnp_dir} \\
        -o vireo_out \\
        -N ${meta.n_donors} \\
        -p ${task.cpus} \\
        --randSeed 42
    """
}

// ─── AMULET ───────────────────────────────────────────────────────────────────
// ATAC doublet detection from fragment overlaps.
// Source: 05_quality_control.sh:336-342, 422-428

process AMULET {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_amulet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC/AMULET" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("amulet_out/MultipletSummary.txt"),         emit: summary
    tuple val(meta), path("amulet_out/MultipletBarcodes.txt"),        emit: barcodes
    path "amulet_out/*"

    script:
    def autosomes = (meta.species == 'human') ? '/opt/AMULET/human_autosomes.txt' : '/opt/AMULET/mouse_autosomes.txt'
    def restriction = (meta.species == 'human') \
        ? '/opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed' \
        : '/opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_mm10.bed'
    """
    fragments=\$(find -L ${outs_dir} -name 'fragments.tsv.gz'   | head -1)
    singlecell=\$(find -L ${outs_dir} -name 'singlecell.csv'    | head -1)

    mkdir -p amulet_out

    AMULET.sh \\
        "\$fragments" \\
        "\$singlecell" \\
        ${autosomes} \\
        ${restriction} \\
        amulet_out \\
        /opt/AMULET/
    """
}

// ─── MGATK2 ───────────────────────────────────────────────────────────────────
// Mitochondrial genotyping for ATAC libraries.
// Source: 05_quality_control.sh:373, 432

process MGATK2 {
    tag "$meta.library_id"
    label 'process_medium'   // overridden to 32c/128GB/96h via withName: 'MGATK2'
    container "${params.container_mgatk}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("mgatk2/"), emit: results

    script:
    def bam       = "${outs_dir}/possorted_bam.bam"
    def barcodes  = "${outs_dir}/filtered_peak_bc_matrix/barcodes.tsv"
    """
    mkdir -p mgatk2

    mgatk2 run \\
        -i  ${bam} \\
        -o  mgatk2 \\
        -b  ${barcodes} \\
        -c  ${task.cpus}
    """
}

// ─── MACS3 ────────────────────────────────────────────────────────────────────
// Custom peak calling on ATAC fragment files.
// Settings follow ENCODE scATAC recommendations:
//   --nomodel --shift -75 --extsize 150 (nucleosome-free region model)
//   --keep-dup all                       (cellranger-atac already deduplicates)
//   --nolambda                           (disable local background; sparse libraries)
// Source: 04_count.sh (cellranger-atac outs), ENCODE ATAC pipeline recommendations

process MACS3 {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_macs3}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("peaks/"), emit: peaks

    script:
    def gsize = (meta.species == 'human') ? 'hs' : 'mm'
    """
    fragments=\$(find -L ${outs_dir} -name 'fragments.tsv.gz' | head -1)
    if [ -z "\$fragments" ]; then
        echo "ERROR: fragments.tsv.gz not found in ${outs_dir}" >&2
        exit 1
    fi

    # Stream fragments to BED format directly using process substitution (prevents writing large temp files to disk)
    mkdir -p peaks
    macs3 callpeak \\
        -t <(zcat "\$fragments" | grep -v '^#' | cut -f1-3) \\
        -f BED \\
        -n ${meta.library_id} \\
        -g ${gsize} \\
        --nomodel \\
        --shift -75 \\
        --extsize 150 \\
        --keep-dup all \\
        --nolambda \\
        -q ${params.macs3_qvalue} \\
        --outdir peaks/
    """
}

// ─── VIRAL_DETECT ─────────────────────────────────────────────────────────────
// Detect viral transcripts using simpleaf (piscem + alevin-fry).
// BAM→FASTQ via bamtofastq v1.4.1; piscem index built from RVDB-nt C-RVDBv31.0.
// bamtofastq writes: {outdir}/{libid}_{n}_{n}_{flowcell}/bamtofastq_S1_L001_R{1,2}_001.fastq.gz
// ALEVIN_FRY_HOME set to PWD/.alevin_fry_home per-task for isolation.

process VIRAL_DETECT {
    tag "$meta.library_id"
    label 'process_high'
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(whitelist), val(simpleaf_chemistry)
    path  viral_index
    path  viral_t2g
    path  bamtofastq_bin

    output:
    tuple val(meta), path("viral/"), emit: counts
    path "versions.yml",                                emit: versions

    script:
    """
    chmod +x ${bamtofastq_bin}
    ./${bamtofastq_bin} \\
        --nthreads ${task.cpus} \\
        --relaxed \\
        ${bam} \\
        fastqs/

    r1=\$(find fastqs -name '*_R1_*.fastq.gz' | sort | paste -sd',')
    r2=\$(find fastqs -name '*_R2_*.fastq.gz' | sort | paste -sd',')

    export ALEVIN_FRY_HOME=\${PWD}/.alevin_fry_home
    mkdir -p "\${ALEVIN_FRY_HOME}"
    simpleaf set-paths

    simpleaf quant \\
        --reads1        "\${r1}" \\
        --reads2        "\${r2}" \\
        --threads       ${task.cpus} \\
        --index         ${viral_index} \\
        --chemistry     ${simpleaf_chemistry} \\
        --t2g-map       ${viral_t2g} \\
        --resolution    cr-like \\
        --unfiltered-pl ${whitelist} \\
        --output        viral

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf --version | sed 's/simpleaf //')
        bamtofastq: \$(./bamtofastq_linux --help 2>&1 | head -1 | sed 's/bamtofastq v//')
    END_VERSIONS
    """
}

// ─── SIMPLEAF_VELOCITY ────────────────────────────────────────────────────────
// Spliced/unspliced quantification for RNA velocity.
// Re-quantifies GEX FASTQs using simpleaf USA mode (spliceu reference).
// Runs after cellbender; uses cellbender barcodes as the permitted list.
// Not run for Flex (probe-based chemistry, no intronic signal).
// R import: fishpond::loadFry("velocity/af_quant", outputFormat = "velocity")

process SIMPLEAF_VELOCITY {
    tag "$meta.library_id"
    label 'process_high'
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy'

    input:
    tuple val(meta), val(gex_fastq_dirs), val(simpleaf_chemistry), path(barcodes)
    path spliceu_index

    output:
    tuple val(meta), path("velocity/"), emit: counts
    path "versions.yml",                emit: versions

    script:
    """
    r1=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R1_*.fastq.gz' 2>/dev/null | sort | paste -sd',')
    r2=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R2_*.fastq.gz' 2>/dev/null | sort | paste -sd',')

    if [ -z "\$r1" ] || [ -z "\$r2" ]; then
        echo "ERROR: no FASTQs found for ${meta.id} in: ${gex_fastq_dirs}" >&2
        exit 1
    fi

    # Cellbender outputs barcodes with -1 suffix (CellRanger format); strip for simpleaf
    sed 's/-1\$//' ${barcodes} > barcodes_clean.txt

    export ALEVIN_FRY_HOME=\${PWD}/.alevin_fry_home
    mkdir -p "\${ALEVIN_FRY_HOME}"
    simpleaf set-paths

    simpleaf quant \\
        --reads1        "\${r1}" \\
        --reads2        "\${r2}" \\
        --threads       ${task.cpus} \\
        --index         ${spliceu_index} \\
        --chemistry     ${simpleaf_chemistry} \\
        --t2g-map       ${spliceu_index}/t2g_3col.tsv \\
        --resolution    cr-like \\
        --unfiltered-pl barcodes_clean.txt \\
        --output        velocity

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf --version | sed 's/simpleaf //')
    END_VERSIONS
    """
}

// ─── SCRUBLET ────────────────────────────────────────────────────────────────
// GEX Doublet detection.
// Source: nf-core/scdownstream doublet detection step

process SCRUBLET {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_scrublet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/scrublet" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(cellbender_h5)

    output:
    tuple val(meta), path("scrublet_out/doublets.csv"), emit: doublets
    path "scrublet_out/*"

    script:
    """
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    mkdir -p scrublet_out

    python << 'PYEOF'
import os
import scanpy as sc
import pandas as pd

adata = sc.read_10x_h5('${cellbender_h5}')
sc.pp.scrublet(adata, expected_doublet_rate=0.08)

df = pd.DataFrame({
    'doublet_score': adata.obs['doublet_score'],
    'is_gex_doublet': adata.obs['predicted_doublet'].astype(bool)
}, index=adata.obs_names)
df.index.name = 'barcode'
df.to_csv('scrublet_out/doublets.csv')
PYEOF
    """
}

