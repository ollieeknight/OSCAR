// ─── CELLRANGER_ATAC ──────────────────────────────────────────────────────────
// Handles: DOGMA-ATAC, Multiome-ATAC, standalone ATAC, ASAP-ATAC.
// Output dir: {library_id}_ATAC  (mirrors current OSCAR naming)
// Source: 04_count.sh:78-167

process CELLRANGER_ATAC {
    tag "$meta.library_id"
    label 'process_high'   // overridden to 64c/128GB/96h via withName
    container "${params.container_cellranger_atac}"
    publishDir { "${params.outdir}/${params.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("${meta.library_id}_ATAC/outs"), emit: outs
    path "versions.yml",                                    emit: versions

    script:
    def reference = meta.species == 'human' ? params.ref_human : params.ref_mouse
    // DOGMA ATAC requires --chemistry ARC-v1; standalone ATAC does not
    def extra_args = (meta.assay == 'DOGMA') ? '--chemistry ARC-v1' : ''
    """
    cellranger-atac count \\
        --id        "${meta.library_id}_ATAC" \\
        --reference "${reference}" \\
        --fastqs    "${fastq_dir}" \\
        --sample    "${meta.id}" \\
        --localcores \$(nproc) \\
        --localmem  ${task.memory.toGiga()} \\
        ${extra_args}

    # Clean up cellranger-atac temp dirs
    rm -rf "${meta.library_id}_ATAC/SC_ATAC_COUNTER_CS" "${meta.library_id}_ATAC/_"*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-atac: \$(cellranger-atac --version 2>&1 | head -1 | sed 's/.* //')
    END_VERSIONS
    """
}
