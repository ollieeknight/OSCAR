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

process SIMPLEAF_VELOCITY {
    tag "$meta.library_id"
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy',
               saveAs: { fn -> fn.startsWith('velocity/af_map/') ? null : fn }

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
