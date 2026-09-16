process CELLSNP_LITE {
    tag "$meta.library_id ($mode)"
    container "${params.container_cellsnp}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${mode == 'atac' ? "${meta.library_id}_ATAC" : meta.library_id}/vireo" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)
    val(mode)

    output:
    tuple val(meta), path("cellsnp_${meta.library_id}/"), emit: vcf
    path "cellsnp_${meta.library_id}/*"

    script:
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

process VIREO {
    tag "$meta.library_id"
    container "${params.container_vireo}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${mode == 'atac' ? "${meta.library_id}_ATAC" : meta.library_id}/vireo" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(cellsnp_dir)
    val(mode)

    output:
    tuple val(meta), path("vireo_out/donor_ids.tsv"), emit: donor_ids
    path "vireo_out/*"

    script:
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
