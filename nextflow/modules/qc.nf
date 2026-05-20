// ─── CELLBENDER ───────────────────────────────────────────────────────────────
// Ambient RNA removal. GPU process.
// Source: 05_quality_control.sh:91-139

process CELLBENDER {
    tag "$meta.library_id"
    label 'process_gpu'
    container "${params.container_cellbender}"
    publishDir "${params.outdir}/${meta.library_id}/cellbender", mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("output.h5"),               emit: h5
    tuple val(meta), path("output_cell_barcodes.csv"), emit: barcodes
    path "versions.yml",                               emit: versions

    script:
    """
    feature_matrix=\$(find ${outs_dir} -name 'raw_feature_bc_matrix.h5' | head -1)
    if [ -z "\$feature_matrix" ]; then
        echo "ERROR: raw_feature_bc_matrix.h5 not found in ${outs_dir}" >&2
        exit 1
    fi

    cellbender remove-background \\
        --cuda \\
        --input  "\$feature_matrix" \\
        --output output.h5 \\
        --cpu-threads \$(nproc) \\
        --checkpoint-mins 10000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(cellbender --version 2>&1 | sed 's/CellBender v//')
    END_VERSIONS
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
    publishDir "${params.outdir}/${out_dir}/cellsnp", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)
    val(mode)       // 'gex' or 'atac'

    output:
    tuple val(meta), path("cellsnp_${meta.library_id}/"), emit: vcf
    path "versions.yml",                                   emit: versions

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
        -p  \$(nproc) \\
        ${umi_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp-lite: \$(cellsnp-lite --version 2>&1 | sed 's/cellsnp-lite //')
    END_VERSIONS
    """
}

// ─── VIREO ────────────────────────────────────────────────────────────────────
// Probabilistic donor demultiplexing.
// Source: 05_quality_control.sh:210-214, 279-283, 365-369

process VIREO {
    tag "$meta.library_id"
    label 'process_medium'   // overridden to 32c/64GB/96h via withName: 'CELLSNP_LITE|VIREO'
    container "${params.container_vireo}"
    publishDir "${params.outdir}/${out_dir}/vireo", mode: 'copy'

    input:
    tuple val(meta), path(cellsnp_dir)
    val(mode)   // 'gex' or 'atac'

    output:
    tuple val(meta), path("donor_ids.tsv"),       emit: donor_ids
    tuple val(meta), path("variant_doublets.tsv"), emit: doublets
    path "versions.yml",                           emit: versions

    script:
    def out_dir = (mode == 'atac') ? "${meta.library_id}_ATAC" : meta.library_id
    """
    vireo \\
        -c ${cellsnp_dir} \\
        -o . \\
        -N ${meta.n_donors} \\
        -p \$(nproc)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vireo: \$(vireo --version 2>&1 | head -1 | sed 's/vireo //')
    END_VERSIONS
    """
}

// ─── AMULET ───────────────────────────────────────────────────────────────────
// ATAC doublet detection from fragment overlaps.
// Source: 05_quality_control.sh:336-342, 422-428

process AMULET {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_amulet}"
    publishDir "${params.outdir}/${meta.library_id}_ATAC/AMULET", mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("MultipletSummary.txt"),         emit: summary
    tuple val(meta), path("MultipletBarcodes.txt"),        emit: barcodes
    path "versions.yml",                                   emit: versions

    script:
    """
    fragments=\$(find ${outs_dir} -name 'fragments.tsv.gz'   | head -1)
    singlecell=\$(find ${outs_dir} -name 'singlecell.csv'    | head -1)

    AMULET.sh \\
        "\$fragments" \\
        "\$singlecell" \\
        /opt/AMULET/human_autosomes.txt \\
        /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed \\
        . \\
        /opt/AMULET/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulet: \$(AMULET.sh --version 2>&1 || echo 'v1.0')
    END_VERSIONS
    """
}

// ─── MGATK2 ───────────────────────────────────────────────────────────────────
// Mitochondrial genotyping for ATAC libraries.
// Source: 05_quality_control.sh:373, 432

process MGATK2 {
    tag "$meta.library_id"
    label 'process_medium'   // overridden to 32c/128GB/96h via withName: 'MGATK2'
    container "${params.container_mgatk}"
    publishDir "${params.outdir}/${meta.library_id}_ATAC/mgatk2", mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("mgatk2_out/"), emit: results
    path "versions.yml",                  emit: versions

    script:
    def bam       = "${outs_dir}/possorted_bam.bam"
    def barcodes  = "${outs_dir}/filtered_peak_bc_matrix/barcodes.tsv"
    """
    mkdir -p mgatk2_out

    mgatk2 run \\
        -i  ${bam} \\
        -o  mgatk2_out \\
        -b  ${barcodes} \\
        -c  \$(nproc)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgatk2: \$(mgatk2 --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
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
    label 'process_low'
    container "${params.container_macs3}"
    publishDir "${params.outdir}/${meta.library_id}_ATAC/peaks", mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("peaks/"), emit: peaks
    path "versions.yml",              emit: versions

    script:
    def gsize = (meta.species == 'human') ? 'hs' : 'mm'
    """
    fragments=\$(find ${outs_dir} -name 'fragments.tsv.gz' | head -1)
    if [ -z "\$fragments" ]; then
        echo "ERROR: fragments.tsv.gz not found in ${outs_dir}" >&2
        exit 1
    fi

    # Convert fragment file to BED (strip header, keep chr/start/end only)
    zcat "\$fragments" | grep -v '^#' | cut -f1-3 > fragments_tmp.bed

    mkdir -p peaks
    macs3 callpeak \\
        -t fragments_tmp.bed \\
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

    rm fragments_tmp.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version 2>&1 | sed 's/macs3 //')
    END_VERSIONS
    """
}

// ─── SCRUBLET ─────────────────────────────────────────────────────────────────
// Computational doublet detection for GEX/CITE/DOGMA libraries.
// Runs on the raw (unfiltered) feature-barcode matrix from cellranger multi.
// GEX features are extracted before doublet simulation so ADT/HTO rows in
// CITE/DOGMA libraries do not confound the gene-expression PCA.
// Expected doublet rate: params.scrublet_doublet_rate (default 0.06 / ~6%).

process SCRUBLET {
    tag "$meta.library_id"
    label 'process_low'
    container "${params.container_scrublet}"
    publishDir "${params.outdir}/${meta.library_id}/scrublet", mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("${meta.library_id}_scrublet_scores.csv"), emit: scores
    path "versions.yml",                                              emit: versions

    script:
    """
    raw_matrix=\$(find ${outs_dir} -type d -name 'raw_feature_bc_matrix' | sort | head -1)
    if [ -z "\$raw_matrix" ]; then
        echo "ERROR: raw_feature_bc_matrix not found in ${outs_dir}" >&2
        exit 1
    fi

    python3 << PYEOF
import scrublet as scr
import scipy.io, pandas as pd, numpy as np, os

matrix_dir = "\$raw_matrix"
matrix   = scipy.io.mmread(os.path.join(matrix_dir, 'matrix.mtx.gz')).T.tocsc()
features = pd.read_csv(os.path.join(matrix_dir, 'features.tsv.gz'), sep='\\t',
                       header=None, names=['id', 'name', 'type'])
barcodes = pd.read_csv(os.path.join(matrix_dir, 'barcodes.tsv.gz'), sep='\\t',
                       header=None, names=['barcode'])

# Keep Gene Expression features only — removes ADT/HTO rows in CITE/DOGMA libs
gex_mask = features['type'] == 'Gene Expression'
matrix   = matrix[:, gex_mask.values]

scrub = scr.Scrublet(matrix, expected_doublet_rate=${params.scrublet_doublet_rate})
scores, doublets = scrub.scrub_doublets(
    min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
)

pd.DataFrame({
    'barcode':           barcodes['barcode'].values,
    'doublet_score':     scores.tolist(),
    'predicted_doublet': doublets.tolist()
}).to_csv('${meta.library_id}_scrublet_scores.csv', index=False)

n_total    = len(barcodes)
n_doublets = int(doublets.sum())
pct        = round(float(doublets.mean()) * 100, 1)
print("Scrublet: " + str(n_doublets) + " doublets of " + str(n_total) + " barcodes (" + str(pct) + "%)")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scrublet: \$(python3 -c 'import scrublet; print(scrublet.__version__)' 2>&1)
    END_VERSIONS
    """
}
