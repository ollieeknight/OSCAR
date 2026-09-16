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
    stage_atac_fastqs.py \\
        --sample-id  ${meta.id} \\
        --min-reads  ${min_reads}

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
