include { VALIDATE_FASTQ; FASTP } from '../modules/fastq_qc'

workflow FASTQ_QC {
    take:
        ch_fastqs   // [meta, fastq_dir_str, [fq_files]]

    main:
        // Every FASTQ gets a pigz -t integrity check, and that check alone
        // gates counting. fastp runs as an independent reporting branch on the
        // same files, so its (slower) per-file QC no longer stands between
        // demux and cellranger. The cost is decompressing R-reads twice: once
        // in pigz, once in fastp. They run concurrently, so it costs CPU, not
        // wall time.
        //
        // fastp still only sees R-reads (R1/R2/R3) above 1MB: index reads
        // carry no useful quality signal, and below 1MB the report is noise.
        def fastp_covers = { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }

        ch_fastqs
            .transpose(by: 2)
            .map { meta, fq_dir, fastq -> [meta, fq_dir, fastq, fastq.name] }
            .set { ch_to_validate }

        VALIDATE_FASTQ(ch_to_validate)

        ch_fastqs
            .flatMap { meta, fq_dir, fq_files ->
                def files = fq_files instanceof List ? fq_files : [fq_files]
                files
                    .findAll { f -> fastp_covers.call(f) }
                    .collect { f -> [meta.run_name, fq_dir, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_fastp_input }

        FASTP(ch_fastp_input)

        // Reassemble per-library lists once validation has passed. This channel
        // gates counting, so it stays a real barrier: a file is only re-emitted
        // after the pigz task that verified it completed.
        VALIDATE_FASTQ.out.fastq
            // Grouped by (meta.id, fastq_dir), so one row per library per
            // flowcell. Grouping on meta.id alone merged every flowcell into
            // one row and then discarded all but the first fastq_dir, which
            // hid the other flowcells from velocity quant (it reads the dir
            // string, not the file list). COUNT_GEX regroups these rows by
            // library_id, so cellranger still sees every flowcell's FASTQs.
            .map { meta, fq_dir, fastq -> [[meta.id, fq_dir], meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { key, metas, _fq_dirs, fastqs -> [metas[0], key[1], fastqs] }
            .set { ch_validated_fastqs }

        FASTP.out.report
            .transpose(by: 2)   // one row per report file; the {json,html} glob emits a list
            .groupTuple(by: [0, 1])
            .set { ch_fastp_reports }

    emit:
        fastqs        = ch_validated_fastqs   // [meta, fastq_dir_string, [fastq_files]]
        fastp_reports = ch_fastp_reports      // [run_name, fastq_dir, [report_files]]
}
