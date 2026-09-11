include { GENERATE_SAMPLESHEET  } from '../modules/demux'
include { BCLCONVERT            } from '../modules/demux'
include { CLEAN_FASTQ_DIR       } from '../modules/demux'
include { FASTQ_QC              } from './fastq_qc'

workflow DEMUX {
    take:
        ch_meta_bcl  // [meta, bcl_dir] pairs — each meta pre-bound to its BCL dir

    main:
        // Wipe stale published fastq.gz from any earlier invocation before this
        // run's BCLCONVERT tasks publish into the same directory (publishDir
        // only adds/overwrites — it never deletes files a previous run left
        // behind under a different bcl-convert-assigned S-number).
        ch_meta_bcl
            .map { meta, bcl_dir -> [bcl_dir.name, bcl_dir.parent.toString()] }
            .unique()
            .set { ch_unique_bcl_dirs }

        CLEAN_FASTQ_DIR(ch_unique_bcl_dirs)

        // Group by (assay_chemistry, bcl_dir) → one samplesheet per group.
        // Deliberately NOT split by modality or index_type: bcl-convert supports
        // mixed SI(8bp)/DI(10bp) rows in one [BCLConvert_Data] table (SI rows
        // just leave Index2 blank — see is_dual/data_rows below), and the
        // OverrideCycles mask table in lib/indexes.nf never actually varies by
        // modality for a given chemistry+index_type. Merging GEX+ADT+HTO into
        // one bcl-convert pass per lane cuts flowcell-tile reads 3x → 1x.
        ch_meta_bcl
            .map { meta, bcl_dir ->
                def key = "${meta.assay}_${meta.chemistry}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs ->
                // Materialise ArrayBag → ArrayList before passing to any process
                def ml      = []
                metas.each { m -> ml << m }
                def bcl_dir = bcl_dirs[0]

                // Defensive: index length must agree within each index kind (SI vs
                // DI). Mixing SI and DI *between* kinds in one group is expected
                // (that's the whole point); two different kit lengths within the
                // same kind would be a real samplesheet error.
                [false, true].each { is_dual_flag ->
                    def subset = ml.findAll { m -> (m.index_seqs?.is_dual ?: false) == is_dual_flag }
                    if (subset.size() > 1) {
                        def len0 = subset[0].index_seqs?.rows[0]?.i7?.length() ?: 0
                        if (subset.any { m -> (m.index_seqs?.rows[0]?.i7?.length() ?: 0) != len0 })
                            error "Demux group ${key} has mixed ${is_dual_flag ? 'DI' : 'SI'} index lengths: " +
                                  subset.collect { m -> "${m.id}=${m.index_seqs?.rows[0]?.i7?.length() ?: 0}" }.join(', ')
                    }
                }

                // Pre-build samplesheet data section — avoids ArrayBag ops inside process script
                def is_dual     = ml.any { m -> m.index_seqs.is_dual }
                def data_header = is_dual ? 'Sample_ID,Index,Index2' : 'Sample_ID,Index'
                def data_rows   = ml.collectMany { m ->
                    m.index_seqs.rows.collect { row ->
                        is_dual ? "${m.id},${row.i7},${row.get('i5', '')}" : "${m.id},${row.i7}"
                    }
                }.join('\n')

                [key, ml, bcl_dir, bcl_dir.parent.toString(), is_dual, data_header, data_rows]
            }
            .set { ch_demux_input }

        GENERATE_SAMPLESHEET(ch_demux_input)

        // Detect which lanes have cbcl data (Groovy filesystem read — same pattern as detect_sequencer).
        // flatMap emits one channel item per present lane → one BCLCONVERT job per lane.
        // This avoids --no-lane-splitting memory buffering on high-output flowcells.
        // Each item is gated on CLEAN_FASTQ_DIR for its bcl_dir via combine(by:0),
        // so no BCLCONVERT task can publish into a fastq dir that hasn't been wiped.
        GENERATE_SAMPLESHEET.out.samplesheet
            .flatMap { demux_key, metas, bcl_dir, bcl_parent, samplesheet ->
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
                present_lanes.collect { lane -> [bcl_dir.name, demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane] }
            }
            .combine(CLEAN_FASTQ_DIR.out.done, by: 0)
            .map { bcl_name, demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane, cleaned ->
                [demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane]
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
        fastqs        = FASTQ_QC.out.fastqs        // [meta, fastq_dir_string, [validated_fastq_files]]
        fastp_reports = FASTQ_QC.out.fastp_reports // [run_name, fastq_dir, [report_files]] (passed to MULTIQC in main.nf)
}
