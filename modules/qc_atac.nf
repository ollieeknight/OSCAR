process AMULET {
    tag "$meta.library_id"
    container "${params.container_amulet}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC/AMULET" }, mode: 'copy',
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(outs_dir)

    output:
    path "amulet/*"

    script:
    def genome = meta.species == 'human' ? 'hg38' : 'mm10'
    """
    mkdir -p amulet

    AMULET.sh \\
        ${outs_dir}/fragments.tsv.gz \\
        ${outs_dir}/singlecell.csv \\
        /opt/AMULET/${meta.species}_autosomes.txt \\
        /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_${genome}.bed \\
        amulet \\
        /opt/AMULET/
    """
}

process MGATK2 {
    tag "$meta.library_id"
    container "${params.container_mgatk}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    path "mgatk2/"

    script:
    """
    mkdir -p mgatk2

    mgatk2 run \\
        -i ${outs_dir}/possorted_bam.bam \\
        -o mgatk2 \\
        -b ${outs_dir}/filtered_peak_bc_matrix/barcodes.tsv \\
        -c ${task.cpus}
    """
}

process MACS3 {
    tag "$meta.library_id"
    container "${params.container_macs3}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC" }, mode: 'copy'

    input:
    tuple val(meta), path(outs_dir)

    output:
    path "peaks/"

    script:
    """
    macs3 callpeak \\
        -t <(zcat ${outs_dir}/fragments.tsv.gz | grep -v '^#' | cut -f1-3) \\
        -f BED \\
        -n ${meta.library_id} \\
        -g ${meta.species == 'human' ? 'hs' : 'mm'} \\
        --nomodel \\
        --shift -75 \\
        --extsize 150 \\
        --keep-dup all \\
        --nolambda \\
        -q ${params.macs3_qvalue} \\
        --outdir peaks
    """
}
