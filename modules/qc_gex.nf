process CELLBENDER {
    tag "$meta.library_id"
    container "${params.container_cellbender}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/cellbender" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("cellbender/output_filtered.h5"),       emit: h5
    tuple val(meta), path("cellbender/output_cell_barcodes.csv"), emit: barcodes
    path "cellbender/*"

    script:
    """
    mkdir -p cellbender

    cellbender remove-background \\
        --cuda \\
        --input           ${outs_dir}/multi/count/raw_feature_bc_matrix.h5 \\
        --output          cellbender/output.h5 \\
        --cpu-threads     ${task.cpus} \\
        --checkpoint-mins 10000
    """
}

process SCRUBLET {
    tag "$meta.library_id"
    container "${params.container_scrublet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/scrublet" }, mode: 'copy'

    input:
    tuple val(meta), path(cellbender_h5)

    output:
    path "doublets.csv"

    script:
    """
    export MPLCONFIGDIR=\$PWD/tmp/mpl NUMBA_CACHE_DIR=\$PWD/tmp/numba
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}

    call_doublets.py \\
        --h5           ${cellbender_h5} \\
        --min-counts   ${params.scrublet_min_counts} \\
        --doublet-rate ${params.scrublet_doublet_rate} \\
        --out          doublets.csv
    """
}
