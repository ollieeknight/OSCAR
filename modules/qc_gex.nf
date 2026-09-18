process CELLBENDER {
    tag "$meta.library_id"
    container "${params.container_cellbender}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/cellbender" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("cellbender_out/output_filtered.h5"),       emit: h5
    tuple val(meta), path("cellbender_out/output_cell_barcodes.csv"), emit: barcodes
    path "cellbender_out/*"

    script:
    """
    feature_matrix=\$(find -L ${outs_dir} -name 'raw_feature_bc_matrix.h5' | sort | head -1)
    if [ -z "\$feature_matrix" ]; then
        echo "ERROR: raw_feature_bc_matrix.h5 not found in ${outs_dir}" >&2
        exit 1
    fi

    mkdir -p cellbender_out

    cellbender remove-background \\
        --cuda \\
        --input  "\$feature_matrix" \\
        --output cellbender_out/output.h5 \\
        --cpu-threads ${task.cpus} \\
        --checkpoint-mins 10000
    """
}

process SCRUBLET {
    tag "$meta.library_id"
    container "${params.container_scrublet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}/scrublet" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(cellbender_h5)

    output:
    tuple val(meta), path("scrublet_out/doublets.csv"), emit: doublets
    path "scrublet_out/*"

    script:
    """
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NUMBA_NUM_THREADS=${task.cpus}

    mkdir -p scrublet_out

    python << 'PYEOF'
import os
import scanpy as sc
import pandas as pd

adata = sc.read_10x_h5('${cellbender_h5}')
sc.pp.scrublet(adata, expected_doublet_rate=0.08)

df = pd.DataFrame({
    'doublet_score': adata.obs['doublet_score'],
    'is_gex_doublet': adata.obs['predicted_doublet'].astype(bool)
}, index=adata.obs_names)
df.index.name = 'barcode'
df.to_csv('scrublet_out/doublets.csv')
PYEOF
    """
}

