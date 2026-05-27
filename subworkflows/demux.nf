include { GENERATE_SAMPLESHEET } from '../modules/demux'
include { BCLCONVERT }        from '../modules/demux'
include { FALCO }               from '../modules/demux'
include { VALIDATE_FASTQ }      from '../modules/demux'

workflow DEMUX {
    take:
        ch_meta_bcl  // [meta, bcl_dir] pairs — each meta pre-bound to its BCL dir

    main:
        // Group by (assay_indextype_chemistry_modality_indexlength, bcl_dir) → one BCL Convert job per group
        // Index length determines OverrideCycles: 8bp (TruSeq) vs 10bp (10x)
        ch_meta_bcl
            .map { meta, bcl_dir ->
                def index_len = meta.index_seqs?.rows[0]?.i7?.length() ?: 0
                def key = "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs ->
                // Materialise ArrayBag → ArrayList before passing to any process
                def ml      = []
                metas.each { m -> ml << m }
                def bcl_dir = bcl_dirs[0]

                // Validate index-length homogeneity (guaranteed by key, but defensive)
                def index_len = ml[0].index_seqs?.rows[0]?.i7?.length() ?: 10
                if (ml.any { m -> (m.index_seqs?.rows[0]?.i7?.length() ?: 10) != index_len })
                    error "Demux group ${key} has mixed index lengths: " +
                          ml.collect { m -> "${m.id}=${m.index_seqs?.rows[0]?.i7?.length() ?: 10}" }.join(', ')

                // Pre-build samplesheet data section — avoids ArrayBag ops inside process script
                def is_dual     = ml.any { m -> m.index_seqs.is_dual }
                def data_header = is_dual ? 'Sample_ID,Index,Index2' : 'Sample_ID,Index'
                def data_rows   = ml.collectMany { m ->
                    m.index_seqs.rows.collect { row ->
                        is_dual ? "${m.id},${row.i7},${row.get('i5', '')}" : "${m.id},${row.i7}"
                    }
                }.join('\n')

                [key, ml, bcl_dir, is_dual, data_header, data_rows]
            }
            .set { ch_demux_input }

        GENERATE_SAMPLESHEET(ch_demux_input)
        BCLCONVERT(GENERATE_SAMPLESHEET.out.samplesheet)

        // Explode output back to individual metas by matching meta.id to FASTQ filename.
        // meta.id == BCL Convert Sample_ID == FASTQ filename prefix.
        // When the same library appears across multiple BCL dirs (multi-run), both
        // sets of FASTQs end up as separate [meta, fastqs] items; groupTuple in main.nf
        // merges them before counting.
        // Derive the published fastq dir path from bcl_dir.name (e.g. R463_bcl → R463_fastq).
        // Pass the directory string rather than staged files — prevents filename collisions
        // when the same library is sequenced on multiple flowcells (identical _S1_ naming).
        BCLCONVERT.out.fastqs
            .flatMap { metas, bcl_name, fq_files ->
                def fqs     = fq_files instanceof List ? fq_files : [fq_files]
                def run     = bcl_name.replaceAll(/_bcl.*$/, '')
                def fq_dir  = fqs[0].parent.toAbsolutePath().toString()
                metas.collectMany { meta ->
                    def matched = fqs.findAll { f -> f.name.contains(meta.id) }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        // Transpose FASTQ file lists to validate each individual file as a separate task
        ch_fastqs
            .transpose(by: 2)
            .map { meta, fq_dir, fastq -> [meta, fq_dir, fastq, fastq.name] }
            .set { ch_to_validate }

        VALIDATE_FASTQ(ch_to_validate)

        // Group validated files back into lists per library meta
        VALIDATE_FASTQ.out.fastq
            .map { meta, fq_dir, fastq -> [meta.id, meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { id, metas, fq_dirs, fastqs ->
                [metas[0], fq_dirs.unique(false), fastqs]
            }
            .set { ch_validated_fastqs }

        // Run Falco per R-read FASTQ (R1/R2/R3 only; I1/I2 index reads skipped)
        // Thread run_name (derived from bcl_dir) so each report lands in the correct
        // {run}_fastq/falco/ directory, not a shared params.run_name dir.
        BCLCONVERT.out.fastqs
            .flatMap { metas, bcl_name, fq_files ->
                def run = bcl_name.replaceAll(/_bcl.*$/, '')
                (fq_files instanceof List ? fq_files : [fq_files])
                    .findAll { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }
                    .collect { f -> [run, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_falco_input }

        FALCO(ch_falco_input)

        FALCO.out.report
            .groupTuple(by: 0)
            .set { ch_falco_reports }

    emit:
        fastqs        = ch_validated_fastqs  // [meta, fastq_dir_string, [validated_fastq_files]]
        falco_reports = ch_falco_reports     // collected falco dirs (passed to MULTIQC in main.nf)
        versions      = BCLCONVERT.out.versions.mix(VALIDATE_FASTQ.out.versions)
}
