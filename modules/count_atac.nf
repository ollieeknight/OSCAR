process CELLRANGER_ATAC {
    tag "$meta.library_id"
    container "${params.container_cellranger_atac}"
    publishDir { "${params.outdir}/${meta.run_name}_outs" }, mode: 'copy'

    input:
    tuple val(meta), path(atac_fastqs, stageAs: 'fastqs/atac/run_???/*')

    output:
    tuple val(meta), path("${meta.library_id}_ATAC/outs"), emit: outs

    script:
    def reference = meta.species == 'human' ? params.ref_human : params.ref_mouse
    def chemistry = meta.assay == 'DOGMA' ? '--chemistry ARC-v1' : ''
    """
    stage_fastqs.py atac

    cellranger-atac count \\
        --id         ${meta.library_id}_ATAC \\
        --reference  ${reference} \\
        --fastqs     staged/atac \\
        --sample     ${meta.id} \\
        --localcores ${task.cpus} \\
        --localmem   ${task.memory.toGiga()} ${chemistry}
    """
}
