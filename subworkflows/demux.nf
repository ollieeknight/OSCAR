include { GENERATE_SAMPLESHEET; BCLCONVERT; CLEAN_FASTQ_DIR; DEMUX_QC } from '../modules/demux'
include { FASTP } from '../modules/fastq_qc'

workflow DEMUX {
    take:
        ch_meta_bcl

    main:
        ch_meta_bcl
            .map { _meta, bcl_dir -> [bcl_dir.name, bcl_dir.parent.toString()] }
            .unique()
            .set { ch_unique_bcl_dirs }

        CLEAN_FASTQ_DIR(ch_unique_bcl_dirs)

        ch_meta_bcl
            .map { meta, bcl_dir ->

                def mergeable = meta.assay in ['CITE', 'GEX']
                def key = mergeable
                    ? "${meta.assay}_${meta.chemistry}_${bcl_dir.name}"
                    : "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs ->
                def meta_list = metas.toSorted { m -> m.id }
                def bcl_dir   = bcl_dirs.toSorted { d -> d.toString() }[0]

                def is_dual = meta_list.any { m -> m.index_seqs.is_dual }

                def sample_specs = meta_list.collectMany { m ->
                    m.index_seqs.rows.collect { row ->
                        [
                            id:         m.id,
                            i7:         row.i7,
                            i5:         row.get('i5', ''),
                            index_len:  row.i7.length(),
                            is_dual:    m.index_seqs.is_dual,
                            assay:      m.assay,
                            chemistry:  m.chemistry,
                            index_type: m.index_type,
                            modality:   m.modality,
                        ]
                    }
                }

                def seen = [:]
                sample_specs.each { sp ->
                    def bc = "${sp.i7}+${sp.i5}"
                    if (seen.containsKey(bc) && seen[bc] != sp.id)
                        error "Demux group ${key}: index ${bc} is used by both " +
                              "'${seen[bc]}' and '${sp.id}'. Two libraries in one lane " +
                              "cannot share a barcode. Check metadata.csv."
                    seen[bc] = sp.id
                }

                [key, meta_list, bcl_dir, bcl_dir.parent.toString(), is_dual, sample_specs]
            }
            .set { ch_demux_input }

        // metas stay out of demux inputs so editing one library's metadata does not re-demultiplex its lane group.
        GENERATE_SAMPLESHEET(ch_demux_input.map { key, _metas, bcl_dir, bcl_parent, is_dual, specs ->
            [key, bcl_dir, bcl_parent, is_dual, specs]
        })

        GENERATE_SAMPLESHEET.out.samplesheet
            .flatMap { demux_key, bcl_dir, bcl_parent, samplesheet ->
                def base_calls = new File("${bcl_dir}/Data/Intensities/BaseCalls")
                def present_lanes = base_calls.listFiles()
                    ?.findAll { d -> d.isDirectory() && d.name =~ /^L\d+$/ }
                    ?.findAll { lane_dir ->
                        new File("${lane_dir}/C1.1").listFiles()?.any { f -> f.name.endsWith('.cbcl') }
                    }
                    ?.collect { d -> d.name.replaceAll(/^L0*/, '').toInteger() }
                    ?.sort()
                if (!present_lanes)
                    error "No lanes with cbcl data found in ${bcl_dir}/Data/Intensities/BaseCalls/"
                present_lanes.collect { lane -> [bcl_dir.name, demux_key, bcl_dir, bcl_parent, samplesheet, lane] }
            }
            .combine(CLEAN_FASTQ_DIR.out.done, by: 0)
            .map { _bcl_name, demux_key, bcl_dir, bcl_parent, samplesheet, lane, _cleaned ->
                [demux_key, bcl_dir, bcl_parent, samplesheet, lane]
            }
            .set { ch_bclconvert_input }

        BCLCONVERT(ch_bclconvert_input)

        BCLCONVERT.out.reports
            .map { _demux_key, bcl_name, bcl_parent, _lane, stats, unknown ->
                def run = bcl_name.replaceAll(/_bcl.*$/, '')
                [run, "${bcl_parent}/${run}_fastq", stats, unknown]
            }
            .groupTuple(by: [0, 1])
            .map { run, fq_dir, stats, unknown ->
                [run, fq_dir,
                 // Every file is named Demultiplex_Stats.csv; sort on the full path.
                 stats.toSorted { s -> s.toString() },
                 unknown.toSorted { u -> u.toString() }]
            }
            .set { ch_demux_qc }

        DEMUX_QC(ch_demux_qc, channel.value(file("${projectDir}/assets/indexes")))

        BCLCONVERT.out.fastqs
            .groupTuple(by: [0, 1])
            .join(ch_demux_input.map { key, metas, _bcl_dir, _bcl_parent, _is_dual, _specs -> [key, metas] })
            .flatMap { demux_key, bcl_name, bcl_parents, fq_file_lists, metas ->
                def fqs     = fq_file_lists.flatten().toSorted { f -> f.name }
                def run     = bcl_name.replaceAll(/_bcl.*$/, '')
                def fq_dir  = "${bcl_parents[0]}/${run}_fastq".toString()

                def valid_fqs = fqs.findAll { f -> f.size() >= 30 }
                def empty_fqs = fqs.findAll { f -> f.size() < 30 }
                if (empty_fqs) log.warn "WARNING: Dropping ${empty_fqs.size()} empty FASTQ(s) from ${demux_key}: ${empty_fqs*.name.join(', ')}"

                metas.collectMany { meta ->
                    def id_re = java.util.regex.Pattern.compile(
                        "^" + java.util.regex.Pattern.quote(meta.id) + "_S\\d+_")
                    def matched = valid_fqs.findAll { f -> id_re.matcher(f.name).find() }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        // Index reads and near-empty files are not worth a fastp report.
        FASTP(ch_fastqs.flatMap { meta, fq_dir, fqs ->
            fqs.findAll { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }
               .collect { f -> [meta.run_name, fq_dir, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
        })

    emit:
        ch_fastqs
}
