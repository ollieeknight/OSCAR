include { VALIDATE_FASTQ; FASTP } from '../modules/fastq_qc'

workflow FASTQ_QC {
    take:
        ch_fastqs   // [meta, fastq_dir_str, [fq_files]]

    main:
        // Integrity checking and fastp both decompress the whole file, so each
        // FASTQ goes to one of them, never both. fastp doubles as the gzip
        // check: on a truncated stream its reader hits
        // error_exit("igzip: unexpected eof"), which is exit(-1), failing the
        // task. (fastp issue #410, the exit-0 bug, is a malformed record on
        // --stdin rather than a truncated file.)
        //
        // Two kinds of file fastp never sees, which pigz checks instead:
        //   - I1/I2 index reads (no _R*_ in the name)
        //   - R-reads at or below the 1MB fastp floor
        def fastp_covers = { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }

        ch_fastqs
            .transpose(by: 2)
            .branch { _meta, _fq_dir, fastq ->
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
        // Keyed on (fastq_dir, basename): bcl-convert's naming has no flowcell
        // component, so two flowcells of the same library+lane produce
        // identical basenames. `join` is 1:1, so a basename-only key silently
        // drops every flowcell after the first — halving a library's reads with
        // no error. failOnDuplicate turns any future key collision into a
        // crash rather than lost data.
        ch_split.via_fastp
            .map { meta, fq_dir, fastq ->
                [ [fq_dir, fastq.name.replaceAll(/\.fastq\.gz$/, '').toString()], meta, fq_dir, fastq ]
            }
            .join(FASTP.out.checked.map { d, n, ok -> [[d, n], ok] }, by: 0, failOnDuplicate: true)
            .map { _key, meta, fq_dir, fastq, _ok -> [meta, fq_dir, fastq] }
            .set { ch_fastp_checked }

        VALIDATE_FASTQ.out.fastq
            .mix(ch_fastp_checked)
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
