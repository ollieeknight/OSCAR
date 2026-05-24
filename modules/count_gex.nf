// ─── MULTI_CONFIG ──────────────────────────────────────────────────────────────
// Generates the cellranger multi CSV from channel metadata and validates that
// each (fastq_dir, meta.id) pair has at least one FASTQ before writing the
// [libraries] rows. Separated from CELLRANGER_MULTI so the heavyweight compute
// process has no Groovy logic and cannot hit ArrayBag serialisation errors.
process MULTI_CONFIG {
    tag "$library_id"
    label 'process_low'

    input:
    tuple val(library_id), val(metas), val(fastq_dirs), path(adt_csv)

    output:
    tuple val(library_id), val(metas), path("multi_config.csv"), emit: config
    path "versions.yml",                                          emit: versions

    script:
    def metas_list = []
    metas.each { m -> if (m != null) metas_list << m }
    if (!metas_list) error "MULTI_CONFIG: metas empty for library ${library_id}"
    def meta = metas_list[0]

    def dirs_list = []
    fastq_dirs.each { d -> dirs_list << d }

    def is_human             = meta.species == 'human'
    def ref_gex              = is_human ? params.ref_human : params.ref_mouse
    def ref_vdj              = is_human ? params.ref_vdj_human : params.ref_vdj_mouse
    def is_dogma_or_multiome = meta.assay in ['DOGMA', 'Multiome']
    def is_flex              = meta.assay == 'Flex'
    def has_vdj              = metas_list.any { it.modality in ['VDJ-T', 'VDJ-B'] }
    def has_adt              = metas_list.any { it.modality in ['ADT', 'HTO'] }
    def create_bam           = (is_human && meta.n_donors > 1) ? 'true' : 'false'

    // Build CSV lines as a plain list — joining with \n gives exact blank-line control
    def csv_lines = []
    csv_lines << '[gene-expression]'
    csv_lines << "reference,${ref_gex}"
    csv_lines << "create-bam,${create_bam}"
    if (is_dogma_or_multiome)            csv_lines << 'chemistry,ARC-v1'
    else if (is_flex && meta.chemistry)  csv_lines << "chemistry,${meta.chemistry}"

    if (has_vdj) {
        csv_lines << ''
        csv_lines << '[vdj]'
        csv_lines << "reference,${ref_vdj}"
    }

    if (has_adt && adt_csv.name != 'NO_FILE') {
        csv_lines << ''
        csv_lines << '[feature]'
        csv_lines << "reference,${adt_csv.toAbsolutePath()}"
    }

    csv_lines << ''
    csv_lines << '[libraries]'
    csv_lines << 'fastq_id,fastqs,feature_types'

    def config_header = csv_lines.join('\n')

    // Emit one library row per (dir, modality) pair only when FASTQs are present.
    // Uses POSIX find -print -quit so the check is fast and container-portable.
    def lib_check_lines = []
    dirs_list.each { dir ->
        metas_list.each { m ->
            def ft = m.modality == 'GEX'            ? 'Gene Expression'      :
                     m.modality in ['ADT', 'HTO']   ? 'Antibody Capture'     :
                     m.modality == 'VDJ-T'          ? 'VDJ-T'                :
                     m.modality == 'VDJ-B'          ? 'VDJ-B'                :
                     m.modality == 'CRISPR'         ? 'CRISPR Guide Capture' : 'Gene Expression'
            lib_check_lines << \
                "find \"${dir}\" -maxdepth 2 -name \"${m.id}*.fastq.gz\" " +
                "-print -quit 2>/dev/null | grep -q . && " +
                "echo \"${m.id},${dir},${ft}\" >> multi_config.csv || true"
        }
    }
    def lib_check_script = lib_check_lines.join('\n')

    """
    cat > multi_config.csv << 'MULTI_CONFIG_EOF'
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
