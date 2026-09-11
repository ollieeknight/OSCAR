include { GENERATE_SAMPLESHEET; BCLCONVERT; CLEAN_FASTQ_DIR; DEMUX_QC } from '../modules/demux'
include { FASTQ_QC } from './fastq_qc'

workflow DEMUX {
    take:
        ch_meta_bcl  // [meta, bcl_dir] pairs — each meta pre-bound to its BCL dir

    main:
        // publishDir only adds and overwrites, so a re-run under a changed
        // samplesheet leaves the old S-numbered fastq alongside the new ones.
        // Wipe the target dir first.
        ch_meta_bcl
            .map { _meta, bcl_dir -> [bcl_dir.name, bcl_dir.parent.toString()] }
            .unique()
            .set { ch_unique_bcl_dirs }

        CLEAN_FASTQ_DIR(ch_unique_bcl_dirs)

        // Group by (assay, chemistry, bcl_dir), so one bcl-convert call per lane
        // covers every modality and index type of a flowcell. A CITE-seq
        // flowcell takes 1 pass over the tiles instead of 3.
        //
        // Mixed 8bp SI and 10bp DI libraries share a call via the per-sample
        // OverrideCycles column in [BCLConvert_Data] (BCL Convert 4.1.5+) —
        // see GENERATE_SAMPLESHEET. Assay and chemistry stay in the key because
        // they fix the Y-read structure, which that column cannot vary.
        ch_meta_bcl
            .map { meta, bcl_dir ->
                // Only CITE/GEX merge across modality. Multiome/DOGMA/ASAP ATAC
                // gives read 3 a different role (Y24 cell barcode where GEX
                // masks index2) and emits a different FASTQ set per sample, so
                // those stay one group per modality.
                def mergeable = meta.assay in ['CITE', 'GEX']
                def key = mergeable
                    ? "${meta.assay}_${meta.chemistry}_${bcl_dir.name}"
                    : "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs ->
                // groupTuple emits an ArrayBag; processes need a plain List.
                def meta_list = metas.collect()
                def bcl_dir   = bcl_dirs[0]

                // Index2 column is present when ANY member is dual-indexed;
                // single-index rows leave it blank and mask i5 in their own
                // OverrideCycles entry.
                def is_dual = meta_list.any { m -> m.index_seqs.is_dual }

                // One spec per samplesheet row. A single-index 10x kit code
                // expands to 4 sequences, so rows outnumber metas.
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

                // Duplicate index within one bcl-convert call is a hard error:
                // bcl-convert cannot assign a read between two samples sharing
                // a barcode, and merging groups is what makes this reachable.
                def seen = [:]
                sample_specs.each { sp ->
                    def bc = "${sp.i7}+${sp.i5}"
                    if (seen.containsKey(bc) && seen[bc] != sp.id)
                        error "Demux group ${key}: index ${bc} is used by both " +
                              "'${seen[bc]}' and '${sp.id}'. Two libraries in one lane " +
                              "cannot share a barcode — check metadata.csv."
                    seen[bc] = sp.id
                }

                [key, meta_list, bcl_dir, bcl_dir.parent.toString(), is_dual, sample_specs]
            }
            .set { ch_demux_input }

        GENERATE_SAMPLESHEET(ch_demux_input)

        // One BCLCONVERT job per lane that actually holds cbcl data, which
        // avoids the --no-lane-splitting memory buffering on high-output
        // flowcells. combine(by:0) gates each job on CLEAN_FASTQ_DIR for its
        // own bcl_dir, so nothing publishes into a dir that is not wiped yet.
        GENERATE_SAMPLESHEET.out.samplesheet
            .flatMap { demux_key, metas, bcl_dir, bcl_parent, samplesheet ->
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
                present_lanes.collect { lane -> [bcl_dir.name, demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane] }
            }
            .combine(CLEAN_FASTQ_DIR.out.done, by: 0)
            .map { _bcl_name, demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane, _cleaned ->
                [demux_key, metas, bcl_dir, bcl_parent, samplesheet, lane]
            }
            .set { ch_bclconvert_input }

        BCLCONVERT(ch_bclconvert_input)

        // Post-demux QC, one job per run. Each group and lane writes its own
        // Reports/, so group by run and concatenate inside the process. Judging
        // a group alone would make a whole missing group look like a small pool.
        BCLCONVERT.out.reports
            .map { _demux_key, _metas, bcl_name, bcl_parent, _lane, stats, unknown ->
                def run = bcl_name.replaceAll(/_bcl.*$/, '')
                [run, "${bcl_parent}/${run}_fastq", stats, unknown]
            }
            .groupTuple(by: [0, 1])
            .set { ch_demux_qc }

        DEMUX_QC(ch_demux_qc, channel.value(file("${projectDir}/assets/indexes")))

        // Merge FASTQs from all lanes of the same demux group, then explode back to individual metas.
        // groupTuple(by: [0,1]) groups by (demux_key, bcl_dir_name) — same group across lanes.
        BCLCONVERT.out.fastqs
            .groupTuple(by: [0, 2])
            .flatMap { demux_key, metas_per_lane, _bcl_name, fq_file_lists ->
                // metas identical across lanes — take first; flatten all lane FASTQs into one list
                def metas   = metas_per_lane[0]
                def fqs     = fq_file_lists.flatten()
                def fq_dir  = fqs[0].parent.toAbsolutePath().toString()

                // Branch: warn on empty FASTQs (<30 bytes), drop from downstream
                def valid_fqs = fqs.findAll { f -> f.size() >= 30 }
                def empty_fqs = fqs.findAll { f -> f.size() < 30 }
                if (empty_fqs) log.warn "WARNING: Dropping ${empty_fqs.size()} empty FASTQ(s) from ${demux_key}: ${empty_fqs*.name.join(', ')}"

                metas.collectMany { meta ->
                    // Anchor on bcl-convert's naming (<Sample_ID>_S<n>_L<lane>_<read>_001)
                    // rather than a bare substring test: merging every modality of a
                    // flowcell into one group puts far more libraries in the same
                    // directory, so an id that is a prefix of another id would
                    // otherwise steal its FASTQ.
                    def id_re = java.util.regex.Pattern.compile(
                        "^" + java.util.regex.Pattern.quote(meta.id) + "_S\\d+_")
                    def matched = valid_fqs.findAll { f -> id_re.matcher(f.name).find() }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        FASTQ_QC(ch_fastqs)

    emit:
        fastqs        = FASTQ_QC.out.fastqs        // [meta, fastq_dir_string, [validated_fastq_files]]
        fastp_reports = FASTQ_QC.out.fastp_reports // [run_name, fastq_dir, [report_files]] (passed to MULTIQC in main.nf)
        demux_summary = DEMUX_QC.out.summary       // [run_name, fastq_dir, summary_csv]
        demux_mqc     = DEMUX_QC.out.mqc           // [run_name, fastq_dir, custom-content CSV for MultiQC]
}
