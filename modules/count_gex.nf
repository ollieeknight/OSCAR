process CELLRANGER_MULTI {
    tag "$library_id"
    container "${params.container_cellranger}"
    publishDir { "${params.outdir}/${run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(run_name), val(config_header), path(adt_csv), val(flex_samples),
          path(gex_fastqs,    stageAs: 'fastqs/gex/run_???/*'),
          path(adt_fastqs,    stageAs: 'fastqs/adt/run_???/*'),
          path(hto_fastqs,    stageAs: 'fastqs/hto/run_???/*'),
          path(vdj_t_fastqs,  stageAs: 'fastqs/vdj_t/run_???/*'),
          path(vdj_b_fastqs,  stageAs: 'fastqs/vdj_b/run_???/*'),
          path(crispr_fastqs, stageAs: 'fastqs/crispr/run_???/*')

    output:
    tuple val(library_id), path("${library_id}/outs"), emit: outs

    script:
    """
    cat > config_header.txt << 'HEADER_EOF'
${config_header}
HEADER_EOF

    cat > flex_samples.txt << 'SAMPLES_EOF'
${flex_samples}
SAMPLES_EOF

    stage_fastqs.py multi

    cellranger multi \\
        --id         ${library_id} \\
        --csv        multi_config.csv \\
        --localcores ${task.cpus} \\
        --localmem   ${task.memory.toGiga()}
    """
}
