process CELLRANGER_MULTI {
    tag "$library_id"
    container "${params.container_cellranger}"
    publishDir { "${params.outdir}/${metas[0].run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas), val(config_header), path(adt_csv),
          val(flex_samples_content),
          path(gex_fastqs,    stageAs: "fastqs/gex/run_???/*"),
          path(adt_fastqs,    stageAs: "fastqs/adt/run_???/*"),
          path(hto_fastqs,    stageAs: "fastqs/hto/run_???/*"),
          path(vdj_t_fastqs,  stageAs: "fastqs/vdj_t/run_???/*"),
          path(vdj_b_fastqs,  stageAs: "fastqs/vdj_b/run_???/*"),
          path(crispr_fastqs, stageAs: "fastqs/crispr/run_???/*")

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs

    script:
    def min_reads = 10000
    """
    cat > config_header.txt << 'OSCAR_CR_HEADER_EOF'
${config_header}
OSCAR_CR_HEADER_EOF

    cat > flex_samples.txt << 'OSCAR_FLEX_SAMPLES_EOF'
${flex_samples_content}
OSCAR_FLEX_SAMPLES_EOF

    stage_multi_fastqs.py \\
        --library-id ${library_id} \\
        --min-reads  ${min_reads}

    cellranger multi \\
        --id        "${library_id}" \\
        --csv       multi_config.csv \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}
    """
}
