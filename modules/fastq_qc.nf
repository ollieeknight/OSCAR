// ─── FASTQ quality control ───────────────────────────────────────────────────
// Per-file QC and integrity checking, plus the run-level MultiQC report.

// ─── FASTP ───────────────────────────────────────────────────────────────────
// Per-file QC, report only. One job per R-read FASTQ (R1/R2/R3); index reads
// (I1/I2) skipped. No -o/--out1, so fastp writes reports and no filtered reads.

process FASTP {
    tag "$fastq_name"
    container "${params.container_fastp}"
    // Published beside the source flowcell's FASTQs (fastq_dir), not under
    // params.outdir — otherwise every --extra_bcl_dirs run drops a stray
    // {run}_fastq/fastp tree into the primary run's directory.
    publishDir { "${fastq_dir}/fastp" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), val(fastq_name), path(fastq)

    output:
    tuple val(run_name), val(fastq_dir), path("${run_name}_${fastq_name}.{json,html}"), emit: report
    // Completion signal keyed by the FASTQ this task verified. FASTQ_QC uses it
    // to release the file downstream, so a fastp failure blocks its library.
    // Keyed on (fastq_dir, fastq_name), not fastq_name alone: bcl-convert puts
    // no flowcell in the filename, so the same library+lane on two flowcells
    // yields identical basenames and the release join stops being 1:1.
    tuple val(fastq_dir), val(fastq_name), val(true), emit: checked

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

// ─── MULTIQC ─────────────────────────────────────────────────────────────────

process MULTIQC {
    container "${params.container_multiqc}"
    // Beside the source flowcell's FASTQs, matching FASTP.
    publishDir { "${fastq_dir}/multiqc" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), path(reports)

    output:
    path "multiqc_report.html",      emit: report
    path "multiqc_report_data/",     emit: data

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    multiqc ${config} --force --filename multiqc_report -o . .
    """
}

// ─── VALIDATE_FASTQ ──────────────────────────────────────────────────────────
// Lightweight validation step that runs gzip -t on each individual fastq file.
// Fully distributed across Slurm nodes and benefits from Nextflow caching.

process VALIDATE_FASTQ {
    tag "$fastq_name"
    container "${params.container_pigz}"

    input:
    tuple val(meta), val(fastq_dir), path(fastq), val(fastq_name)

    output:
    tuple val(meta), val(fastq_dir), path(fastq), emit: fastq

    script:
    """
    pigz -t -f -p ${task.cpus} ${fastq}
    """
}
