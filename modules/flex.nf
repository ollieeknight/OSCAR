include { get_flex_barcode_file; get_flex_whitelist_file } from '../lib/chemistry'

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

process FLEX_BARCODE_EXTRACT {
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "probe_barcodes.txt", emit: barcodes

    script:
    def bc_file = get_flex_barcode_file(chemistry)
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/translation/${bc_file} probe_barcodes.txt
    """
}

process FLEX_WHITELIST_EXTRACT {
    container "${params.container_cellranger}"

    input:
    val(chemistry)

    output:
    path "cb_whitelist.txt.gz", emit: whitelist

    script:
    def wl_file = get_flex_whitelist_file(chemistry)
    """
    cp /opt/cellranger-10.0.0/lib/python/cellranger/barcodes/${wl_file} cb_whitelist.txt.gz
    """
}

process FLEX_SAMPLE_PREPARE {
    container "${params.container_python}"

    input:
    path samples_file
    path probe_barcodes_ref

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
    publishDir { "${params.outdir}/${metas[0].run_name}_outs" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas),
          path(probe_tsv_cyto),
          path(cyto_probe_barcodes),
          path(cb_whitelist),
          val(cyto_preset),
          path(gex_fastqs, stageAs: "fastqs/gex/run_???/*")

    output:
    tuple val(library_id), val(metas), path("${library_id}_cyto"), emit: counts

    script:
    def min_reads = 10000
    """
    stage_cyto_fastqs.py --min-reads ${min_reads}

    FASTQ_PAIRS=\$(cat fastq_pairs.txt)

    PROBES_ARG=""
    if [ "${cyto_probe_barcodes}" != "NO_FILE" ]; then
        PROBES_ARG="-p ${cyto_probe_barcodes}"
    fi

    cyto workflow gex \\
        -c ${probe_tsv_cyto} \\
        \${PROBES_ARG} \\
        -w ${cb_whitelist} \\
        --preset ${cyto_preset} \\
        -o ${library_id}_cyto \\
        -F mtx \\
        --no-filter \\
        --memory-limit ${params.flex_cyto_memory_limit} \\
        -T ${task.cpus} \\
        -f \\
        \${FASTQ_PAIRS}
    """
}

process CYTO_RENAME_SAMPLES {
    tag "$library_id"

    input:
    tuple val(library_id), val(metas), path(cyto_out), path(samples_file)

    output:
    tuple val(library_id), val(metas), path(cyto_out), emit: counts

    script:
    """
    cyto_rename_samples.py \\
        --samples-file ${samples_file} \\
        --cyto-out     ${cyto_out}
    """
}
