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
    publishDir { "${params.outdir}/${params.run_name}_outs/${library_id}" }, mode: 'copy'

    input:
    tuple val(library_id), val(metas), val(fastq_dirs), path(adt_csv)

    output:
    tuple val(library_id), val(metas), path("${library_id}/outs"), emit: outs
    path "versions.yml",                                            emit: versions

    script:
    def meta                 = metas[0]   // all metas share species, n_donors, etc.
    def is_human             = meta.species == 'human'
    def ref_gex              = is_human ? params.ref_human : params.ref_mouse
    def ref_vdj              = is_human ? params.ref_vdj_human : params.ref_vdj_mouse
    def is_dogma_or_multiome = meta.assay in ['DOGMA', 'Multiome']
    def is_flex              = meta.assay == 'Flex'

    def has_vdj    = metas.any { it.modality in ['VDJ-T', 'VDJ-B'] }
    def has_adt    = metas.any { it.modality in ['ADT', 'HTO'] }
    def has_crispr = metas.any { it.modality == 'CRISPR' }

    // BAM file only needed for donor demultiplexing (cellsnp-lite/vireo).
    // Skip for mouse (no human SNP VCF) and single-donor runs — saves 20-50GB
    // of storage per library and 20-30% of Cell Ranger runtime.
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

    // fastq_id = meta.id (= BCL Convert Sample_ID = filename prefix)
    // fastqs   = published fastq dir — one row per (modality, flowcell) combination.
    // Listing the same fastq_id twice with different dirs is how cellranger multi
    // handles libraries sequenced across multiple flowcells.
    // Entries where the sample's FASTQs in a given dir are < 10 MB are skipped —
    // this handles modalities absent from a flowcell that still produce near-empty files.
    def dirs_list = (fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]) as ArrayList
    def lib_checks = []
    dirs_list.each { dir ->
        metas.each { m ->
            def ft = (m.modality == 'GEX')           ? 'Gene Expression' :
                     (m.modality in ['ADT', 'HTO'])   ? 'Antibody Capture' :
                     (m.modality == 'VDJ-T')          ? 'VDJ-T' :
                     (m.modality == 'VDJ-B')          ? 'VDJ-B' :
                     (m.modality == 'CRISPR')         ? 'CRISPR Guide Capture' : 'Gene Expression'
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
