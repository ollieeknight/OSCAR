// ─── Supplementary quantification ────────────────────────────────────────────
// Optional re-quantification of GEX data: viral transcript detection and
// spliced/unspliced counts for RNA velocity. Both use simpleaf.

// ─── VIRAL_DETECT ────────────────────────────────────────────────────────────
// Detect viral transcripts using simpleaf (piscem + alevin-fry).
// BAM→FASTQ via bamtofastq v1.4.1; piscem index built from RVDB-nt C-RVDBv31.0.
// bamtofastq writes: {outdir}/{libid}_{n}_{n}_{flowcell}/bamtofastq_S1_L001_R{1,2}_001.fastq.gz
// ALEVIN_FRY_HOME set to PWD/.alevin_fry_home per-task for isolation.

process VIRAL_DETECT {
    tag "$meta.library_id"
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(whitelist), val(simpleaf_chemistry)
    path  viral_index
    path  viral_t2g
    path  bamtofastq_bin

    output:
    tuple val(meta), path("viral/"), emit: counts

    script:
    """
    chmod +x ${bamtofastq_bin}
    ./${bamtofastq_bin} \\
        --nthreads ${task.cpus} \\
        --relaxed \\
        ${bam} \\
        fastqs/

    r1=\$(find fastqs -name '*_R1_*.fastq.gz' | sort | paste -sd',')
    r2=\$(find fastqs -name '*_R2_*.fastq.gz' | sort | paste -sd',')

    export ALEVIN_FRY_HOME=\${PWD}/.alevin_fry_home
    mkdir -p "\${ALEVIN_FRY_HOME}"
    simpleaf set-paths

    simpleaf quant \\
        --reads1        "\${r1}" \\
        --reads2        "\${r2}" \\
        --threads       ${task.cpus} \\
        --index         ${viral_index} \\
        --chemistry     ${simpleaf_chemistry} \\
        --t2g-map       ${viral_t2g} \\
        --resolution    cr-like \\
        --unfiltered-pl ${whitelist} \\
        --output        viral
    """
}

// ─── SIMPLEAF_VELOCITY ───────────────────────────────────────────────────────
// Spliced/unspliced quantification for RNA velocity.
// Re-quantifies GEX FASTQs using simpleaf USA mode (spliceu reference).
// Runs after cellbender; uses cellbender barcodes as the permitted list.
// Not run for Flex (probe-based chemistry, no intronic signal).
// R import: fishpond::loadFry("velocity/af_quant", outputFormat = "velocity")

process SIMPLEAF_VELOCITY {
    tag "$meta.library_id"
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy'

    input:
    tuple val(meta), val(gex_fastq_dirs), val(simpleaf_chemistry), path(barcodes)
    path spliceu_index

    output:
    tuple val(meta), path("velocity/"), emit: counts

    script:
    """
    r1=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R1_*.fastq.gz' 2>/dev/null | sort | paste -sd',')
    r2=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R2_*.fastq.gz' 2>/dev/null | sort | paste -sd',')

    if [ -z "\$r1" ] || [ -z "\$r2" ]; then
        echo "ERROR: no FASTQs found for ${meta.id} in: ${gex_fastq_dirs}" >&2
        exit 1
    fi

    # Cellbender outputs barcodes with -1 suffix (CellRanger format); strip for simpleaf
    sed 's/-1\$//' ${barcodes} > barcodes_clean.txt

    export ALEVIN_FRY_HOME=\${PWD}/.alevin_fry_home
    mkdir -p "\${ALEVIN_FRY_HOME}"
    simpleaf set-paths

    simpleaf quant \\
        --reads1        "\${r1}" \\
        --reads2        "\${r2}" \\
        --threads       ${task.cpus} \\
        --index         ${spliceu_index} \\
        --chemistry     ${simpleaf_chemistry} \\
        --t2g-map       ${spliceu_index}/t2g_3col.tsv \\
        --resolution    cr-like \\
        --unfiltered-pl barcodes_clean.txt \\
        --output        velocity
    """
}
