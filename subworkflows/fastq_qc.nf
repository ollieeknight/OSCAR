include { VALIDATE_FASTQ } from '../modules/demux'
include { FASTP          } from '../modules/demux'

workflow FASTQ_QC {
    take:
        ch_fastqs   // [meta, fastq_dir_str, [fq_files]]

    main:
        // Validate each FASTQ individually as a separate task
        ch_fastqs
            .transpose(by: 2)
            .map { meta, fq_dir, fastq -> [meta, fq_dir, fastq, fastq.name] }
            .set { ch_to_validate }

        VALIDATE_FASTQ(ch_to_validate)

        // Reassemble validated files into per-library lists
        VALIDATE_FASTQ.out.fastq
            .map { meta, fq_dir, fastq -> [meta.id, meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { id, metas, fq_dirs, fastqs ->
                // fq_dirs is one entry per validated file, all identical for a given
                // meta.id — collapse to the single dir string downstream expects.
                [metas[0], fq_dirs.unique(false).first(), fastqs]
            }
            .set { ch_validated_fastqs }

        // fastp QC on R-reads only (R1/R2/R3; I1/I2 skipped)
        ch_fastqs
            .flatMap { meta, fq_dir, fq_files ->
                def files = fq_files instanceof List ? fq_files : [fq_files]
                files
                    .findAll { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }
                    .collect { f -> [meta.run_name, fq_dir, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_fastp_input }

        FASTP(ch_fastp_input)

        FASTP.out.report
            .groupTuple(by: [0, 1])
            .set { ch_fastp_reports }

    emit:
        fastqs        = ch_validated_fastqs   // [meta, fastq_dir_string, [fastq_files]]
        fastp_reports = ch_fastp_reports      // [run_name, fastq_dir, [report_files]]
}
