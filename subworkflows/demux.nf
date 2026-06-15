include { GENERATE_SAMPLESHEET } from '../modules/demux'
include { BCLCONVERT }          from '../modules/demux'
include { FASTQ_QC }            from './fastq_qc'

workflow DEMUX {
    take:
        ch_meta_bcl  // [meta, bcl_dir] pairs — each meta pre-bound to its BCL dir

    main:
        // Group by (assay_indextype_chemistry_modality_indexlength, bcl_dir) → one samplesheet per group
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

        // Detect which lanes have cbcl data (Groovy filesystem read — same pattern as detect_sequencer).
        // flatMap emits one channel item per present lane → one BCLCONVERT job per lane.
        // This avoids --no-lane-splitting memory buffering on high-output flowcells.
        GENERATE_SAMPLESHEET.out.samplesheet
            .flatMap { demux_key, metas, bcl_dir, samplesheet ->
                def base_calls = new File("${bcl_dir}/Data/Intensities/BaseCalls")
                def present_lanes = base_calls.listFiles()
                    ?.findAll { it.isDirectory() && it.name =~ /^L\d+$/ }
                    ?.findAll { lane_dir ->
                        new File("${lane_dir}/C1.1").listFiles()?.any { it.name.endsWith('.cbcl') }
                    }
                    ?.collect { it.name.replaceAll(/^L0*/, '').toInteger() }
                    ?.sort()
                if (!present_lanes)
                    error "No lanes with cbcl data found in ${bcl_dir}/Data/Intensities/BaseCalls/"
                present_lanes.collect { lane -> [demux_key, metas, bcl_dir, samplesheet, lane] }
            }
            .set { ch_bclconvert_input }

        BCLCONVERT(ch_bclconvert_input)

        // Merge FASTQs from all lanes of the same demux group, then explode back to individual metas.
        // groupTuple(by: [0,1]) groups by (demux_key, bcl_dir_name) — same group across lanes.
        BCLCONVERT.out.fastqs
            .groupTuple(by: [0, 2])
            .flatMap { demux_key, metas_per_lane, bcl_name, fq_file_lists ->
                // metas identical across lanes — take first; flatten all lane FASTQs into one list
                def metas   = metas_per_lane[0]
                def fqs     = fq_file_lists.flatten()
                def fq_dir  = fqs[0].parent.toAbsolutePath().toString()

                // Branch: warn on empty FASTQs (<30 bytes), drop from downstream
                def valid_fqs = fqs.findAll { f -> f.size() >= 30 }
                def empty_fqs = fqs.findAll { f -> f.size() < 30 }
                if (empty_fqs) log.warn "WARNING: Dropping ${empty_fqs.size()} empty FASTQ(s) from ${demux_key}: ${empty_fqs*.name.join(', ')}"

                metas.collectMany { meta ->
                    def matched = valid_fqs.findAll { f -> f.name.contains(meta.id) }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        FASTQ_QC(ch_fastqs)

    emit:
        fastqs        = FASTQ_QC.out.fastqs        // [meta, [fastq_dir_strings], [validated_fastq_files]]
        falco_reports = FASTQ_QC.out.falco_reports // [run_name, [report_dirs]] (passed to REPORT in main.nf)
}
