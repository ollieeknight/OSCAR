include { BCL_TO_FASTQ } from '../modules/demux'
include { FALCO }        from '../modules/demux'

workflow DEMUX {
    take:
        ch_meta     // channel of meta maps (one per samplesheet row)
        ch_bcl_dirs // channel of BCL dir paths (1 per flowcell; multiple for multi-run)

    main:
        // Group samplesheet rows by demux key: {assay}_{index_type}_{chemistry}_{modality}
        // Each unique combination → one BCL Convert job per BCL dir
        ch_meta
            .map { meta ->
                def key = "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}"
                [key, meta]
            }
            .groupTuple(by: 0)
            .set { ch_demux_groups }

        // Cross-product: each demux group × each BCL dir.
        // Append BCL dir name to key to prevent work-dir collision in multi-run.
        ch_demux_groups
            .combine(ch_bcl_dirs)
            .map { key, metas, bcl_dir ->
                ["${key}_${bcl_dir.name}", metas, bcl_dir]
            }
            .set { ch_demux_input }

        BCL_TO_FASTQ(ch_demux_input)

        // Explode output back to individual metas by matching meta.id to FASTQ filename.
        // meta.id == BCL Convert Sample_ID == FASTQ filename prefix.
        // When the same library appears across multiple BCL dirs (multi-run), both
        // sets of FASTQs end up as separate [meta, fastqs] items; groupTuple in main.nf
        // merges them before counting.
        BCL_TO_FASTQ.out.fastqs
            .flatMap { metas, fq_files ->
                def fqs = fq_files instanceof List ? fq_files : [fq_files]
                metas.collectMany { meta ->
                    def matched = fqs.findAll { f -> f.name.contains(meta.id) }
                    matched ? [[meta, matched]] : []
                }
            }
            .set { ch_fastqs }

        // Run Falco on all FASTQs per demux job for MultiQC
        BCL_TO_FASTQ.out.fastqs
            .map { metas, fq_files ->
                def key = "${metas[0].assay}_${metas[0].index_type}_${metas[0].chemistry}_${metas[0].modality}"
                [key, fq_files instanceof List ? fq_files : [fq_files]]
            }
            .set { ch_falco_input }

        FALCO(ch_falco_input)

        FALCO.out.report
            .collect()
            .set { ch_falco_reports }

    emit:
        fastqs        = ch_fastqs        // [meta, [fastq_files]]
        falco_reports = ch_falco_reports // collected falco dirs (passed to MULTIQC in main.nf)
        versions      = BCL_TO_FASTQ.out.versions
}
