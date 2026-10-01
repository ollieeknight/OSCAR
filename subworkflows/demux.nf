include { GENERATE_SAMPLESHEET; BCLCONVERT; CLEAN_FASTQ_DIR; DEMUX_QC } from '../modules/demux'
include { FASTP } from '../modules/fastq_qc'
include { fastq_dir_for; run_name_for; sample_fastqs } from '../lib/indexes'

workflow DEMUX {
    take:
        ch_meta_bcl

    main:
        CLEAN_FASTQ_DIR(ch_meta_bcl.map { _meta, bcl_dir -> fastq_dir_for(bcl_dir) }.unique())

        ch_demux_input = ch_meta_bcl
            .map { meta, bcl_dir ->
                def mergeable = meta.assay in ['CITE', 'GEX']
                def key = mergeable
                    ? "${meta.assay}_${meta.chemistry}_${bcl_dir.name}"
                    : "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple()
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
                        error "demux group ${key}: '${seen[bc]}' and '${sp.id}' share index ${bc}"
                    seen[bc] = sp.id
                }

                [key, meta_list, bcl_dir, fastq_dir_for(bcl_dir), is_dual, sample_specs]
            }

        // metas stay out of demux inputs so editing one library's metadata does not re-demultiplex its lane group.
        GENERATE_SAMPLESHEET(ch_demux_input.map { key, _metas, bcl_dir, fastq_dir, is_dual, specs ->
            [key, bcl_dir, fastq_dir, is_dual, specs]
        })

        ch_bclconvert_input = GENERATE_SAMPLESHEET.out.samplesheet
            .flatMap { demux_key, bcl_dir, fastq_dir, samplesheet ->
                def lanes = new File("${bcl_dir}/Data/Intensities/BaseCalls").listFiles()
                    ?.findAll { d -> d.name =~ /^L\d+$/ && new File(d, 'C1.1').listFiles()?.any { f -> f.name.endsWith('.cbcl') } }
                    ?.collect { d -> d.name.replaceAll(/^L0*/, '').toInteger() }
                    ?.sort()
                if (!lanes)
                    error "no lanes with cbcl data found in ${bcl_dir}/Data/Intensities/BaseCalls/"
                lanes.collect { lane -> [fastq_dir, demux_key, bcl_dir, samplesheet, lane] }
            }
            .combine(CLEAN_FASTQ_DIR.out.done, by: 0)
            .map { fastq_dir, demux_key, bcl_dir, samplesheet, lane -> [demux_key, bcl_dir, fastq_dir, samplesheet, lane] }

        BCLCONVERT(ch_bclconvert_input)

        ch_demux_qc = BCLCONVERT.out.reports
            .groupTuple()
            .map { fastq_dir, stats, unknown ->
                [run_name_for(file(fastq_dir)), fastq_dir,
                 // Every file is named Demultiplex_Stats.csv; sort on the full path.
                 stats.toSorted { s -> s.toString() },
                 unknown.toSorted { u -> u.toString() }]
            }

        DEMUX_QC(ch_demux_qc, channel.value(file("${projectDir}/assets/indexes")))

        // Index reads, Undetermined and near-empty files are not worth a fastp report.
        FASTP(BCLCONVERT.out.fastqs.flatMap { _key, fastq_dir, fqs ->
            fqs.findAll { f -> f.name =~ /_R\d+_/ && !f.name.startsWith('Undetermined') && f.size() > 1024 * 1024 }
               .collect { f -> [run_name_for(file(fastq_dir)), fastq_dir, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
        })

    emit:
        BCLCONVERT.out.fastqs
            .groupTuple(by: [0, 1])
            .join(ch_demux_input.map { key, metas, _bcl_dir, _fastq_dir, _is_dual, _specs -> [key, metas] })
            .flatMap { demux_key, fastq_dir, fq_file_lists, metas ->
                def (fqs, empty_fqs) = fq_file_lists.flatten().toSorted { f -> f.name }.split { f -> f.size() >= 30 }
                if (empty_fqs) log.warn "dropping ${empty_fqs.size()} empty FASTQ(s) from ${demux_key}: ${empty_fqs*.name.join(', ')}"

                metas.collectMany { meta ->
                    def matched = sample_fastqs(meta.id, fqs)
                    matched ? [[meta, fastq_dir, matched]] : []
                }
            }
}
