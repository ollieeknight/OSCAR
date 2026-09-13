// ─── ATAC quality control ────────────────────────────────────────────────────
// AMULET doublet detection, mgatk2 mitochondrial genotyping, MACS3 peak calling.

// ─── AMULET ──────────────────────────────────────────────────────────────────
// ATAC doublet detection from fragment overlaps.
// Ported from original/bash/05_quality_control.sh

process AMULET {
    tag "$meta.library_id"
    container "${params.container_amulet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC/AMULET" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("amulet_out/MultipletSummary.txt"),         emit: summary
    tuple val(meta), path("amulet_out/MultipletBarcodes.txt"),        emit: barcodes
    path "amulet_out/*"

    script:
    def autosomes = (meta.species == 'human') ? '/opt/AMULET/human_autosomes.txt' : '/opt/AMULET/mouse_autosomes.txt'
    def restriction = (meta.species == 'human') \
        ? '/opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed' \
        : '/opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_mm10.bed'
    """
    fragments=\$(find -L ${outs_dir} -name 'fragments.tsv.gz'   | head -1)
    singlecell=\$(find -L ${outs_dir} -name 'singlecell.csv'    | head -1)

    mkdir -p amulet_out

    AMULET.sh \\
        "\$fragments" \\
        "\$singlecell" \\
        ${autosomes} \\
        ${restriction} \\
        amulet_out \\
        /opt/AMULET/
    """
}

// ─── MGATK2 ──────────────────────────────────────────────────────────────────
// Mitochondrial genotyping for ATAC libraries.
// Ported from original/bash/05_quality_control.sh

process MGATK2 {
    tag "$meta.library_id"
    container "${params.container_mgatk}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("mgatk2/"), emit: results

    script:
    def bam       = "${outs_dir}/possorted_bam.bam"
    def barcodes  = "${outs_dir}/filtered_peak_bc_matrix/barcodes.tsv"
    """
    mkdir -p mgatk2

    mgatk2 run \\
        -i  ${bam} \\
        -o  mgatk2 \\
        -b  ${barcodes} \\
        -c  ${task.cpus}
    """
}

// ─── MACS3 ───────────────────────────────────────────────────────────────────
// Custom peak calling on ATAC fragment files.
// Settings follow ENCODE scATAC recommendations:
//   --nomodel --shift -75 --extsize 150 (nucleosome-free region model)
//   --keep-dup all                       (cellranger-atac already deduplicates)
//   --nolambda                           (disable local background; sparse libraries)
// Ported from original/bash/04_count.sh; settings per the ENCODE ATAC pipeline

process MACS3 {
    tag "$meta.library_id"
    container "${params.container_macs3}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    tuple val(meta), path("peaks/"), emit: peaks

    script:
    def gsize = (meta.species == 'human') ? 'hs' : 'mm'
    """
    fragments=\$(find -L ${outs_dir} -name 'fragments.tsv.gz' | head -1)
    if [ -z "\$fragments" ]; then
        echo "ERROR: fragments.tsv.gz not found in ${outs_dir}" >&2
        exit 1
    fi

    # Stream fragments to BED format directly using process substitution (prevents writing large temp files to disk)
    mkdir -p peaks
    macs3 callpeak \\
        -t <(zcat "\$fragments" | grep -v '^#' | cut -f1-3) \\
        -f BED \\
        -n ${meta.library_id} \\
        -g ${gsize} \\
        --nomodel \\
        --shift -75 \\
        --extsize 150 \\
        --keep-dup all \\
        --nolambda \\
        -q ${params.macs3_qvalue} \\
        --outdir peaks/
    """
}
