include { VALIDATE_FASTQ } from '../modules/demux'
include { FALCO          } from '../modules/demux'

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
                [metas[0], fq_dirs.unique(false), fastqs]
            }
            .set { ch_validated_fastqs }

        // FALCO QC on R-reads only (R1/R2/R3; I1/I2 skipped)
        ch_fastqs
            .flatMap { meta, fq_dir, fq_files ->
                def files = fq_files instanceof List ? fq_files : [fq_files]
                files
                    .findAll { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }
                    .collect { f -> [meta.run_name, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_falco_input }

        FALCO(ch_falco_input)

        FALCO.out.report
            .groupTuple(by: 0)
            .set { ch_falco_reports }

    emit:
        fastqs        = ch_validated_fastqs   // [meta, [fastq_dir_strings], [fastq_files]]
        falco_reports = ch_falco_reports      // [run_name, [report_dirs]]
}
