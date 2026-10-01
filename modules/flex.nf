process FLEX_PROBE_PREPARE {
    container "${params.container_python}"

    input:
    tuple path(standard_probe_csv), path(custom_probe_csv)

    output:
    path "merged_probes_cellranger.csv", emit: probe_csv_cr
    path "merged_probes_cyto.tsv",       emit: probe_tsv_cyto

    script:
    """
    flex_probe_prepare.py \\
        --standard-probe-csv ${standard_probe_csv} \\
        --custom-probe-csv   ${custom_probe_csv}
    """
}

process FLEX_SAMPLE_PREPARE {
    container "${params.container_python}"

    input:
    path(samples_file)
    path(probe_barcodes_ref)

    output:
    path "cyto_probe_barcodes.txt", emit: cyto_barcodes

    script:
    """
    flex_sample_prepare.py \\
        --probe-barcodes-ref ${probe_barcodes_ref} \\
        --samples-file       ${samples_file}
    """
}

process CYTO_FLEX {
    tag "$library_id"
    container "${params.container_cyto}"

    input:
    tuple val(library_id), val(run_name), path(gex_fastqs, stageAs: 'fastqs/gex/run_???/*'),
          path(probe_tsv), path(probe_barcodes), path(whitelist), val(preset)

    output:
    tuple val(library_id), val(run_name), path("${library_id}_cyto"), emit: counts

    script:
    def probes_arg = probe_barcodes.name == 'NO_FILE' ? '' : "-p ${probe_barcodes}"
    """
    stage_fastqs.py cyto

    cyto workflow gex \\
        -c ${probe_tsv} ${probes_arg} \\
        -w ${whitelist} \\
        --preset ${preset} \\
        -o ${library_id}_cyto \\
        -F mtx \\
        --no-filter \\
        --memory-limit ${params.flex_cyto_memory_limit} \\
        -T ${task.cpus} \\
        -f \\
        \$(cat fastq_pairs.txt)
    """
}

// Renames a copy: moving inside the staged input would edit CYTO_FLEX's cached output.
process CYTO_RENAME_SAMPLES {
    tag "$library_id"
    container "${params.container_python}"
    publishDir { "${params.outdir}/${run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(run_name), path(cyto_out, stageAs: 'raw/*'), path(samples_file)

    output:
    path "${library_id}_cyto"

    script:
    """
    cp -rL ${cyto_out} ${library_id}_cyto
    cyto_rename_samples.py \\
        --samples-file ${samples_file} \\
        --cyto-out     ${library_id}_cyto
    """
}
