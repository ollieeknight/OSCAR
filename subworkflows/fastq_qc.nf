include { VALIDATE_FASTQ; FASTP } from '../modules/fastq_qc'

workflow FASTQ_QC {
    take:
        ch_fastqs

    main:
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

        VALIDATE_FASTQ.out.fastq
            .map { meta, fq_dir, fastq -> [[meta.id, fq_dir], meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { key, metas, _fq_dirs, fastqs ->
                [metas.toSorted { m -> m.toString() }[0], key[1], fastqs.toSorted { f -> f.name }]
            }
            .set { ch_validated_fastqs }

        FASTP.out.report
            .transpose(by: 2)
            .groupTuple(by: [0, 1])
            .set { ch_fastp_reports }

    emit:
        fastqs        = ch_validated_fastqs
        fastp_reports = ch_fastp_reports
}
