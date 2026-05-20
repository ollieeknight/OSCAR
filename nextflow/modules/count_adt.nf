// ─── ASAP-seq ADT counting pipeline ──────────────────────────────────────────
// ASAP-only. Triggered after CELLRANGER_ATAC for ASAP libraries.
// Chain: FEATUREMAP → KALLISTO_INDEX → ASAP_TO_KITE → KALLISTO_BUS
//        → BUSTOOLS_CORRECT → BUSTOOLS_SORT → BUSTOOLS_COUNT
// Source: 04_count.sh:246-342

// ─── FEATUREMAP ───────────────────────────────────────────────────────────────

process FEATUREMAP {
    tag "$meta.library_id"
    label 'process_low'
    container "${params.container_asap}"

    input:
    tuple val(meta), path(adt_csv)

    output:
    tuple val(meta), path("FeaturesMismatch.t2g"), path("FeaturesMismatch.fa"), emit: index_files
    path "versions.yml", emit: versions

    script:
    """
    featuremap ${adt_csv} \\
        --t2g FeaturesMismatch.t2g \\
        --fa  FeaturesMismatch.fa  \\
        --header --quiet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featuremap: \$(featuremap --version 2>&1 || echo 'unknown')
    END_VERSIONS
    """
}

// ─── KALLISTO_INDEX ───────────────────────────────────────────────────────────

process KALLISTO_INDEX {
    tag "$meta.library_id"
    label 'process_low'
    container "${params.container_kallisto}"

    input:
    tuple val(meta), path(t2g), path(fa)

    output:
    tuple val(meta), path(t2g), path("FeaturesMismatch.idx"), emit: index
    path "versions.yml", emit: versions

    script:
    """
    kallisto index -i FeaturesMismatch.idx -k 15 ${fa}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version 2>&1 | head -1 | sed 's/kallisto, version //')
    END_VERSIONS
    """
}

// ─── ASAP_TO_KITE ─────────────────────────────────────────────────────────────
// Converts ATAC-barcode-geometry FASTQs → GEX-barcode-geometry FASTQs.

process ASAP_TO_KITE {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_asap}"

    input:
    tuple val(meta), path(adt_fastqs)

    output:
    tuple val(meta), path("kite_converted/"), emit: converted_fastqs
    path "versions.yml",                      emit: versions

    script:
    def fastq_dirs = adt_fastqs instanceof List \
        ? adt_fastqs.collect { it.parent }.unique().join(',') \
        : adt_fastqs.parent.toString()
    def sample_names = adt_fastqs instanceof List \
        ? adt_fastqs.collect { it.simpleName.replaceAll(/_S[0-9]+.*/, '') }.unique().join(',') \
        : adt_fastqs.simpleName.replaceAll(/_S[0-9]+.*/, '')
    """
    mkdir -p kite_converted

    asap_to_kite \\
        -ff "${fastq_dirs}" \\
        -sp "${sample_names}" \\
        -of kite_converted \\
        -on "${meta.library_id}_ADT" \\
        -c  \$(nproc)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        asap_to_kite: \$(asap_to_kite --version 2>&1 || echo 'unknown')
    END_VERSIONS
    """
}

// ─── KALLISTO_BUS ─────────────────────────────────────────────────────────────

process KALLISTO_BUS {
    tag "$meta.library_id"
    label 'process_high'   // overridden via withName: 'KALLISTO_BUS'
    container "${params.container_kallisto}"

    input:
    tuple val(meta), path(t2g), path(idx), path(converted_dir)

    output:
    tuple val(meta), path(t2g), path("bus_output/"), emit: bus
    path "versions.yml",                              emit: versions

    script:
    """
    mkdir -p bus_output

    kallisto bus \\
        -i ${idx} \\
        -o bus_output \\
        -x 0,0,16:0,16,26:1,0,0 \\
        -t \$(nproc) \\
        ${converted_dir}/${meta.library_id}_ADT/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version 2>&1 | head -1 | sed 's/kallisto, version //')
    END_VERSIONS
    """
}

// ─── BUSTOOLS_CORRECT ─────────────────────────────────────────────────────────

process BUSTOOLS_CORRECT {
    tag "$meta.library_id"
    label 'process_low'
    container "${params.container_bustools}"

    input:
    tuple val(meta), path(t2g), path(bus_dir)
    path(whitelist)

    output:
    tuple val(meta), path(t2g), path(bus_dir), path("output_corrected.bus"), emit: corrected
    path "versions.yml",                                                       emit: versions

    script:
    """
    bustools correct \\
        -w ${whitelist} \\
        ${bus_dir}/output.bus \\
        -o output_corrected.bus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version 2>&1 | head -1 | sed 's/bustools //')
    END_VERSIONS
    """
}

// ─── BUSTOOLS_SORT ────────────────────────────────────────────────────────────

process BUSTOOLS_SORT {
    tag "$meta.library_id"
    label 'process_medium'
    container "${params.container_bustools}"

    input:
    tuple val(meta), path(t2g), path(bus_dir), path(corrected_bus)

    output:
    tuple val(meta), path(t2g), path(bus_dir), path("output_sorted.bus"), emit: sorted
    path "versions.yml",                                                    emit: versions

    script:
    """
    bustools sort \\
        -t \$(nproc) \\
        -o output_sorted.bus \\
        ${corrected_bus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version 2>&1 | head -1 | sed 's/bustools //')
    END_VERSIONS
    """
}

// ─── BUSTOOLS_COUNT ───────────────────────────────────────────────────────────

process BUSTOOLS_COUNT {
    tag "$meta.library_id"
    label 'process_low'
    container "${params.container_bustools}"
    publishDir { "${params.outdir}/${params.run_name}_outs/${meta.library_id}_ATAC/ADT" }, mode: 'copy'

    input:
    tuple val(meta), path(t2g), path(bus_dir), path(sorted_bus)

    output:
    tuple val(meta), path("cells_x_genes*"), emit: counts
    path "versions.yml",                      emit: versions

    script:
    """
    bustools count \\
        -o cells_x_genes \\
        --genecounts \\
        -g ${t2g} \\
        -e ${bus_dir}/matrix.ec \\
        -t ${bus_dir}/transcripts.txt \\
        ${sorted_bus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bustools: \$(bustools version 2>&1 | head -1 | sed 's/bustools //')
    END_VERSIONS
    """
}
