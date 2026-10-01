include { get_override_cycles } from '../lib/override_cycles'

process GENERATE_SAMPLESHEET {
    tag "$demux_key"
    container "${params.container_python}"

    input:
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), val(is_dual), val(sample_specs)

    output:
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), path("SampleSheet.csv"), emit: samplesheet

    script:
    def specs = sample_specs.collect { sp ->
        def oc4 = get_override_cycles(sp, 4)
        [sp.id, sp.i7, sp.i5 ?: '-', oc4, get_override_cycles(sp, 3) ?: oc4].join('\t')
    }.join('\n')
    """
    cat > specs.tsv << 'SPECS_EOF'
${specs}
SPECS_EOF

    make_samplesheet.py \\
        --run-info ${bcl_dir}/RunInfo.xml \\
        --specs    specs.tsv ${is_dual ? '--dual' : ''}
    """
}

// Clears FASTQs from an earlier demux so the published folder holds only this run's files.
process CLEAN_FASTQ_DIR {
    tag "$fastq_dir"
    executor 'local'

    input:
    val(fastq_dir)

    output:
    val(fastq_dir), emit: done

    script:
    """
    mkdir -p "${fastq_dir}"
    rm -f "${fastq_dir}"/*.fastq.gz
    """
}

process BCLCONVERT {
    tag "$demux_key L$lane"
    container "${params.container_bclconvert}"
    publishDir { fastq_dir }, mode: 'copy', pattern: "fastqs/*.fastq.gz", saveAs: { fn -> file(fn).name }
    publishDir { "${fastq_dir}/Reports/${demux_key}_L${lane}" }, mode: 'copy', pattern: "fastqs/Reports/*", saveAs: { fn -> file(fn).name }

    input:
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), path(samplesheet), val(lane)

    output:
    tuple val(demux_key), val(fastq_dir), path("fastqs/*.fastq.gz"), emit: fastqs
    tuple val(fastq_dir), path("fastqs/Reports/Demultiplex_Stats.csv"), path("fastqs/Reports/Top_Unknown_Barcodes.csv"), emit: reports

    script:
    def n_tiles = Math.max(1, (task.cpus / 8).toInteger())
    """
    bcl-convert \\
        --bcl-input-directory           ${bcl_dir} \\
        --output-directory              fastqs \\
        --sample-sheet                  ${samplesheet} \\
        --bcl-only-lane                 ${lane} \\
        --bcl-num-parallel-tiles        ${n_tiles} \\
        --bcl-num-conversion-threads    ${n_tiles} \\
        --bcl-num-compression-threads   ${Math.max(1, (task.cpus / 2).toInteger())} \\
        --bcl-num-decompression-threads ${Math.max(1, (task.cpus / 4).toInteger())}
    """
}

process DEMUX_QC {
    tag "$run_name"
    container "${params.container_python}"
    publishDir { fastq_dir }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), path(stats, stageAs: 'stats/*/'), path(unknown, stageAs: 'unknown/*/')
    path(indexes_dir)

    output:
    path "${run_name}_*.csv"

    script:
    """
    demux_qc.py \\
        --stats         stats/*/*.csv \\
        --unknown       unknown/*/*.csv \\
        --indexes       ${indexes_dir} \\
        --prefix        ${run_name} \\
        --dropout-ratio ${params.demux_qc_dropout_ratio} \\
        --unknown-pct   ${params.demux_qc_unknown_pct}
    """
}
