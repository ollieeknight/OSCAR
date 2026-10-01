process CELLSNP_LITE {
    tag "$meta.library_id ($mode)"
    container "${params.container_cellsnp}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${mode == 'atac' ? "${meta.library_id}_ATAC" : meta.library_id}/vireo" }, mode: 'copy',
               pattern: 'cellsnp/*', saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)
    val(mode)

    output:
    tuple val(meta), path("cellsnp/"), emit: vcf
    path "cellsnp/*"

    script:
    def umi_flag   = (mode == 'atac') ? '--UMItag None' : ''
    """
    mkdir -p cellsnp

    cellsnp-lite \\
        -s  ${bam} \\
        -b  ${barcodes} \\
        -O  cellsnp \\
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
    tuple val(meta), path(cellsnp_dir), val(n_donors)
    val(mode)

    output:
    path "vireo_out/*"

    script:
    """
    mkdir -p vireo_out
    vireo \\
        -c ${cellsnp_dir} \\
        -o vireo_out \\
        -N ${n_donors} \\
        -p ${task.cpus} \\
        --randSeed 42
    """
}
