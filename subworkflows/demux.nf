include { GENERATE_SAMPLESHEET } from '../modules/demux'
include { BCLCONVERT }          from '../modules/demux'
include { FASTQ_QC }            from './fastq_qc'

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
        // Pass the directory string rather than staged files — prevents filename collisions
        // when the same library is sequenced on multiple flowcells (identical _S1_ naming).
        BCLCONVERT.out.fastqs
            .flatMap { metas, bcl_name, fq_files ->
                def fqs    = fq_files instanceof List ? fq_files : [fq_files]
                def fq_dir = fqs[0].parent.toAbsolutePath().toString()
                metas.collectMany { meta ->
                    def matched = fqs.findAll { f -> f.name.contains(meta.id) }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        FASTQ_QC(ch_fastqs)

    emit:
        fastqs        = FASTQ_QC.out.fastqs        // [meta, [fastq_dir_strings], [validated_fastq_files]]
        falco_reports = FASTQ_QC.out.falco_reports // [run_name, [report_dirs]] (passed to REPORT in main.nf)
}
