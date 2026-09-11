include { VALIDATE_FASTQ } from '../modules/demux'
include { FASTP          } from '../modules/demux'

workflow FASTQ_QC {
    take:
        ch_fastqs   // [meta, fastq_dir_str, [fq_files]]

    main:
        // Integrity checking and fastp both have to decompress the whole file,
        // so each FASTQ is routed to exactly ONE of them rather than both.
        //
        // fastp is a real gzip integrity check: on a truncated stream its
        // reader calls error_exit("igzip: unexpected eof"), and error_exit is
        // exit(-1), so the task fails. (The known exit-0 bug, fastp issue #410,
        // is a malformed-record case on --stdin, not truncated file input.)
        // Running `pigz -t` over the same bytes first only re-reads them.
        //
        // Two kinds of file fastp never sees, and which therefore still need an
        // explicit check:
        //   - I1/I2 index reads (no _R*_ in the name)
        //   - R-reads at or below the 1MB fastp floor
        // Everything else is covered by the FASTP task below.
        def fastp_covers = { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }

        ch_fastqs
            .transpose(by: 2)
            .branch { meta, fq_dir, fastq ->
                needs_pigz: !fastp_covers.call(fastq)
                via_fastp:  true
            }
            .set { ch_split }

        ch_split.needs_pigz
            .map { meta, fq_dir, fastq -> [meta, fq_dir, fastq, fastq.name] }
            .set { ch_to_validate }

        VALIDATE_FASTQ(ch_to_validate)

        // fastp QC on R-reads only (R1/R2/R3; I1/I2 skipped)
        ch_fastqs
            .flatMap { meta, fq_dir, fq_files ->
                def files = fq_files instanceof List ? fq_files : [fq_files]
                files
                    .findAll { f -> fastp_covers.call(f) }
                    .collect { f -> [meta.run_name, fq_dir, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_fastp_input }

        FASTP(ch_fastp_input)

        // Reassemble per-library lists once BOTH checks have passed.
        //
        // This channel gates counting, so it must stay a real barrier: a file
        // is only re-emitted after the task that verified it completed. Files
        // checked by fastp are re-emitted keyed on the fastp report for the
        // same file, so a failing fastp task still blocks its library.
        ch_split.via_fastp
            .map { meta, fq_dir, fastq ->
                [ fastq.name.replaceAll(/\.fastq\.gz$/, '').toString(), meta, fq_dir, fastq ]
            }
            .join(FASTP.out.checked, by: 0)
            .map { key, meta, fq_dir, fastq, ok -> [meta, fq_dir, fastq] }
            .set { ch_fastp_checked }

        VALIDATE_FASTQ.out.fastq
            .mix(ch_fastp_checked)
            .map { meta, fq_dir, fastq -> [meta.id, meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { id, metas, fq_dirs, fastqs ->
                // fq_dirs is one entry per validated file, all identical for a given
                // meta.id — collapse to the single dir string downstream expects.
                [metas[0], fq_dirs.unique(false).first(), fastqs]
            }
            .set { ch_validated_fastqs }

        FASTP.out.report
            .transpose(by: 2)   // one row per report file; the {json,html} glob emits a list
            .groupTuple(by: [0, 1])
            .set { ch_fastp_reports }

    emit:
        fastqs        = ch_validated_fastqs   // [meta, fastq_dir_string, [fastq_files]]
        fastp_reports = ch_fastp_reports      // [run_name, fastq_dir, [report_files]]
}
