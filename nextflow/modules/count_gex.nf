// ─── CELLRANGER_MULTI ─────────────────────────────────────────────────────────
// Handles: GEX, CITE, DOGMA-GEX, Multiome-GEX.
// Input: [library_id, metas_list, fastq_files (flat), adt_csv (or NO_FILE)]
// All FASTQ files are staged in the work dir; cellranger multi finds them via
// fastq_id matching (fastq_id = meta.id = BCL Convert Sample_ID = filename prefix).
// Source logic: functions.sh:385-515, 04_count.sh:351-401

process CELLRANGER_MULTI {
    tag "$library_id"
    label 'process_high'   // overridden to 64c/128GB/48h via withName
    container "${params.container_cellranger}"
    publishDir "${params.outdir}/${library_id}", mode: 'copy'

    input:
    tuple val(library_id), val(metas), path(fastq_files), path(adt_csv)

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs
    path "versions.yml",                                            emit: versions

    script:
    def meta                 = metas[0]   // all metas share species, n_donors, etc.
    def is_human             = meta.species == 'human'
    def ref_gex              = is_human ? params.ref_human : params.ref_mouse
    def ref_vdj              = is_human ? params.ref_vdj_human : params.ref_vdj_mouse
    def is_dogma_or_multiome = meta.assay in ['DOGMA', 'Multiome']

    def has_vdj    = metas.any { it.modality in ['VDJ-T', 'VDJ-B'] }
    def has_adt    = metas.any { it.modality in ['ADT', 'HTO'] }
    def has_crispr = metas.any { it.modality == 'CRISPR' }

    def ge_section = """[gene-expression]
reference,${ref_gex}
create-bam,true${is_dogma_or_multiome ? '\nchemistry,ARC-v1' : ''}
"""

    def vdj_section = has_vdj ? """
[vdj]
reference,${ref_vdj}
""" : ""

    // adt_csv is staged in the work dir when provided; name='NO_FILE' means absent
    def feature_section = (has_adt && adt_csv.name != 'NO_FILE') ? """
[feature]
reference,${adt_csv}
""" : ""

    // fastq_id = meta.id (= BCL Convert Sample_ID = filename prefix)
    // fastqs   = '.' (work dir, where all FASTQ files are staged)
    def lib_lines = metas.collect { m ->
        def ft = (m.modality == 'GEX')           ? 'Gene Expression' :
                 (m.modality in ['ADT', 'HTO'])   ? 'Antibody Capture' :
                 (m.modality == 'VDJ-T')          ? 'VDJ-T' :
                 (m.modality == 'VDJ-B')          ? 'VDJ-B' :
                 (m.modality == 'CRISPR')         ? 'CRISPR Guide Capture' : 'Gene Expression'
        "${m.id},.,${ft}"
    }.join('\n')

    """
    cat > multi_config.csv << 'MULTIEOF'
${ge_section}${vdj_section}${feature_section}
[libraries]
fastq_id,fastqs,feature_types
${lib_lines}
MULTIEOF

    cellranger multi \\
        --id        "${library_id}" \\
        --csv       multi_config.csv \\
        --localcores \$(nproc) \\
        --localmem  ${task.memory.toGiga()}

    rm -rf "${library_id}/SC_MULTI_CS" "${library_id}/_"*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | head -1 | sed 's/.* //')
    END_VERSIONS
    """
}
