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
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_counts=${params.scrublet_min_counts})

rate = ${params.expected_doublet_rate}
sc.pp.scrublet(adata, expected_doublet_rate=rate)

# Scrublet's automatic threshold assumes a bimodal simulated-doublet
# histogram. When it is unimodal the auto-pick lands in the tail and calls
# almost nothing, so fall back to calling the expected rate by quantile.
score = adata.obs['doublet_score']
auto_thr = adata.uns['scrublet']['threshold']
if (score > auto_thr).mean() < rate / 2:
    thr = score.quantile(1 - rate)
else:
    thr = auto_thr

df = pd.DataFrame({
    'doublet_score': score,
    'is_gex_doublet': (score > thr).astype(bool)
}, index=adata.obs_names)
df.index.name = 'barcode'
df.loc[df['doublet_score'].isna(), 'is_gex_doublet'] = False
df.to_csv('scrublet_out/doublets.csv')
PYEOF
    """
}

