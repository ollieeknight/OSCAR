// ─── CELLRANGER_MULTI ─────────────────────────────────────────────────────────
process CELLRANGER_MULTI {
    tag "$library_id"
    label 'process_high'
    container "${params.container_cellranger}"
    publishDir { "${params.outdir}/${params.run_name}_outs/${library_id}" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas), val(fastq_dirs), path(adt_csv)

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs
    path "versions.yml",                                            emit: versions

    script:
    def metas_list   = metas.collect { it }
    def meta         = metas_list[0]
    def is_human     = meta.species == 'human'
    def ref_gex      = is_human ? params.ref_human : params.ref_mouse
    def ref_vdj      = is_human ? params.ref_vdj_human : params.ref_vdj_mouse
    def is_dogma_or_multiome = meta.assay in ['DOGMA', 'Multiome']
    def is_flex      = meta.assay == 'Flex'

    def has_vdj    = metas_list.any { it.modality in ['VDJ-T', 'VDJ-B'] }
    def has_adt    = metas_list.any { it.modality in ['ADT', 'HTO'] }
    def has_crispr = metas_list.any { it.modality == 'CRISPR' }

    def create_bam = (is_human && meta.n_donors > 1) ? 'true' : 'false'

    def chem_line   = is_dogma_or_multiome ? '\nchemistry,ARC-v1'
                    : is_flex              ? '\nchemistry,' + meta.chemistry
                    : ''
    def ge_section = """[gene-expression]
reference,${ref_gex}
create-bam,${create_bam}${chem_line}
"""

    def vdj_section = has_vdj ? """
[vdj]
reference,${ref_vdj}
""" : ""

    def feature_section = (has_adt && adt_csv.name != 'NO_FILE') ? """
[feature]
reference,${adt_csv.toAbsolutePath()}
""" : ""

    def dirs_list = (fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]) as ArrayList
    def lib_checks = []
    dirs_list.each { dir ->
        metas_list.each { m ->
            def ft = (m.modality == 'GEX')          ? 'Gene Expression'      :
                     (m.modality in ['ADT', 'HTO']) ? 'Antibody Capture'     :
                     (m.modality == 'VDJ-T')        ? 'VDJ-T'                :
                     (m.modality == 'VDJ-B')        ? 'VDJ-B'                :
                     (m.modality == 'CRISPR')       ? 'CRISPR Guide Capture' : 'Gene Expression'
            lib_checks << """\
if [ \$(find "${dir}" -maxdepth 2 -name "${m.id}*.fastq.gz" -printf '%s\\n' 2>/dev/null | awk '{s+=\$1} END{printf "%.0f\\n", s}') -ge 10485760 ]; then
    echo "${m.id},${dir},${ft}" >> multi_config.csv
fi"""
        }
    }
    def lib_check_script = lib_checks.join('\n')

    """
    cat > multi_config.csv << 'MULTIEOF'
${ge_section}${vdj_section}${feature_section}
[libraries]
fastq_id,fastqs,feature_types
MULTIEOF
${lib_check_script}

    cellranger multi \\
        --id        "${library_id}" \\
        --csv       multi_config.csv \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}

    rm -rf "${library_id}/SC_MULTI_CS" "${library_id}/_"*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | head -1 | sed 's/.* //')
END_VERSIONS
    """
}