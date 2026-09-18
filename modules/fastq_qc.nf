process FASTP {
    tag "$fastq_name"
    container "${params.container_fastp}"
    publishDir { "${fastq_dir}/fastp" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), val(fastq_name), path(fastq)

    output:
    tuple val(run_name), val(fastq_dir), path("${run_name}_${fastq_name}.{json,html}"), emit: report

    script:
    """
    fastp \\
        --in1                       ${fastq} \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --thread                    ${task.cpus} \\
        --json                      ${run_name}_${fastq_name}.json \\
        --html                      ${run_name}_${fastq_name}.html
    """
}

