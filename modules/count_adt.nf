process FEATUREMAP {
    tag "$meta.library_id"
    container "${params.container_asap}"

    input:
    tuple val(meta), path(adt_csv)

    output:
    tuple val(meta), path("features.t2g"), path("features.fa"), emit: features

    script:
    """
    featuremap ${adt_csv} \\
        --t2g features.t2g \\
        --fa  features.fa \\
        --header --quiet
    """
}

// asap_to_kite pairs each -ff folder with one -sp prefix, so FASTQs are linked into sets/<flowcell>/<sample>/:
// flowcells reuse file names, and ADT and HTO are separate samples.
process ASAP_TO_KITE {
    tag "$meta.library_id"
    container "${params.container_asap}"

    input:
    tuple val(meta), path(adt_fastqs, stageAs: 'fastqs/run_???/*')

    output:
    tuple val(meta), path("kite/"), emit: fastqs

    script:
    """
    for fq in fastqs/*/*.fastq.gz; do
        flowcell=\$(head -n1 <(gzip -dc "\$fq") | cut -d: -f3)
        sample=\$(basename "\$fq" | sed -E 's/_S[0-9]+_.*//')
        mkdir -p "sets/\$flowcell/\$sample"
        ln -s "\$PWD/\$fq" "sets/\$flowcell/\$sample/"
    done
    sets=(sets/*/*)
    samples=("\${sets[@]##*/}")

    asap_to_kite \\
        -ff "\$(IFS=,; echo "\${sets[*]}")" \\
        -sp "\$(IFS=,; echo "\${samples[*]}")" \\
        -of kite \\
        -on ADT \\
        -c  ${task.cpus}
    """
}

process KALLISTO_BUS {
    tag "$meta.library_id"
    container "${params.container_kallisto}"

    input:
    tuple val(meta), path(t2g), path(fa), path(kite_dir)

    output:
    tuple val(meta), path(t2g), path("bus/"), emit: bus

    script:
    """
    kallisto index -i features.idx -k 15 ${fa}

    kallisto bus \\
        -i features.idx \\
        -o bus \\
        -x 0,0,16:0,16,26:1,0,0 \\
        -t ${task.cpus} \\
        ${kite_dir}/ADT_R1.fastq.gz \\
        ${kite_dir}/ADT_R2.fastq.gz
    """
}

process BUSTOOLS {
    tag "$meta.library_id"
    container "${params.container_bustools}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}_ATAC/ADT" }, mode: 'copy'

    input:
    tuple val(meta), path(t2g), path(bus_dir)
    path(whitelist)

    output:
    path "cells_x_genes*"

    script:
    """
    bustools correct -w ${whitelist} -p ${bus_dir}/output.bus \\
        | bustools sort -t ${task.cpus} -m ${task.memory.toGiga()}G -o sorted.bus -

    bustools count \\
        -o cells_x_genes \\
        --genecounts \\
        -g ${t2g} \\
        -e ${bus_dir}/matrix.ec \\
        -t ${bus_dir}/transcripts.txt \\
        sorted.bus
    """
}
