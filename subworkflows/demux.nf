include { GENERATE_SAMPLESHEET  } from '../modules/demux'
include { BCLCONVERT            } from '../modules/demux'
include { CLEAN_FASTQ_DIR       } from '../modules/demux'
include { DEMUX_QC              } from '../modules/demux'
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

        // Group by (assay, chemistry, bcl_dir) — one bcl-convert call per lane
        // covering every modality and index type of a flowcell, rather than one
        // per (index_type, modality). For a CITE-seq flowcell that is 1 pass
        // instead of 3, each of which otherwise re-reads the same tiles.
        //
        // Mixing 8bp SI and 10bp DI libraries in one call is done with a
        // per-sample OverrideCycles column in [BCLConvert_Data] (BCL Convert
        // 4.1.5+; 4.5.4 in use here) — see GENERATE_SAMPLESHEET. An earlier
        // attempt to put 8bp indexes under one global 10bp mask failed:
        // bcl-convert matches Index against the declared cycle count exactly
        // and does not prefix-match.
        //
        // index_type and modality stay out of the key; assay and chemistry stay
        // in, because they determine the Y-read structure (Read1/Read2 cycle
        // layout), which is NOT expressible per-sample in a way that mixes
        // safely across assays.
        ch_meta_bcl
            .map { meta, bcl_dir ->
                // Merging is limited to assays whose modalities share one
                // Y-read structure and differ only in index length/type --
                // CITE/GEX (GEX + ADT + HTO), which is the wasteful case.
                //
                // Multiome/DOGMA/ASAP are NOT merged across modality: their
                // ATAC mask assigns a different ROLE to a read (ATAC read 3 is
                // a Y24 cell-barcode read where GEX masks index2), so those
                // modalities emit a different set of FASTQ per sample and are
                // not safely expressible as one per-sample column. They keep
                // the original one-group-per-modality behaviour.
                def mergeable = meta.assay in ['CITE', 'GEX']
                def key = mergeable
                    ? "${meta.assay}_${meta.chemistry}_${bcl_dir.name}"
                    : "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs ->
                // Materialise ArrayBag → ArrayList before passing to any process
                def ml      = []
                metas.each { m -> ml << m }
                def bcl_dir = bcl_dirs[0]

                // Index2 column is present when ANY member is dual-indexed;
                // single-index rows leave it blank and mask i5 in their own
                // OverrideCycles entry.
                def is_dual = ml.any { m -> m.index_seqs.is_dual }

                // One spec per samplesheet row. A single-index 10x kit code
                // expands to 4 sequences, so rows outnumber metas.
                def sample_specs = ml.collectMany { m ->
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

                [key, ml, bcl_dir, bcl_dir.parent.toString(), is_dual, sample_specs]
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

        // Post-demux QC, one job per run. Every demux group and lane of a run
        // contributes its own Reports/, so they are grouped by run here and
        // concatenated inside the process: judging each group in isolation
        // would make a whole missing group look like a normal small pool.
        BCLCONVERT.out.reports
            .map { demux_key, metas, bcl_name, bcl_parent, lane, stats, unknown ->
                def run = bcl_name.replaceAll(/_bcl.*$/, '')
                [run, "${bcl_parent}/${run}_fastq", stats, unknown]
            }
            .groupTuple(by: [0, 1])
            .set { ch_demux_qc }

        DEMUX_QC(ch_demux_qc, Channel.value(file("${projectDir}/assets/indexes")))

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
