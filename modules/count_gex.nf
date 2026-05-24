process MULTI_CONFIG {
    tag "$library_id"
    label 'process_low'

    input:
    tuple val(library_id), val(metas), val(config_header), val(lib_check_script), path(adt_csv)

    output:
    tuple val(library_id), val(metas), path("multi_config.csv"), emit: config
    path "versions.yml",                                          emit: versions

    script:
    """
    cat > multi_config.csv << MULTI_CONFIG_EOF
${config_header}
MULTI_CONFIG_EOF
    ${lib_check_script}

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+')
END_VERSIONS
    """
}

// ─── CELLRANGER_MULTI ─────────────────────────────────────────────────────────
process CELLRANGER_MULTI {
    tag "$library_id"
    label 'process_high'
    container "${params.container_cellranger}"
    publishDir { "${params.outdir}/${params.run_name}_outs/${library_id}" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas), path(multi_config)

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs
    path "versions.yml",                                            emit: versions

    script:
    """
    cellranger multi \\
        --id        "${library_id}" \\
        --csv       "${multi_config}" \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}

    rm -rf "${library_id}/SC_MULTI_CS" "${library_id}/_"*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | head -1 | sed 's/.* //')
END_VERSIONS
    """
}
