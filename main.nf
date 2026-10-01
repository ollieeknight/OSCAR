#!/usr/bin/env nextflow

include { DEMUX }       from './subworkflows/demux'
include { COUNT_GEX }   from './subworkflows/count_gex'
include { COUNT_ATAC }  from './subworkflows/count_atac'
include { COUNT_ADT }   from './subworkflows/count_adt'
include { QC_GEX }      from './subworkflows/qc_gex'
include { QC_ATAC }     from './subworkflows/qc_atac'
include { QUANT_EXTRA } from './subworkflows/quant_extra'

include { load_si_indexes; detect_sequencer; run_name_for; sample_fastqs } from './lib/indexes'
include { parse_samplesheet; find_meta_conflicts } from './lib/samplesheet'
include { resolve_input_paths; resolve_run_from; resolve_extras; preflight } from './lib/preflight'

def counts_with_multi(meta) {
    (meta.modality in ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] && meta.assay != 'ASAP') || meta.assay == 'Flex'
}

// [meta, fastq_dir, fastqs] -> [library_id, entries], entries sorted by flowcell then id so -resume sees the same order.
def by_library(ch) {
    ch.map { meta, fastq_dir, fqs -> [meta.library_id, [meta: meta, dir: fastq_dir, files: fqs.toSorted { f -> f.name }]] }
        .groupTuple()
        .map { lid, entries -> [lid, entries.toSorted { a, b -> a.dir <=> b.dir ?: a.meta.id <=> b.meta.id }] }
}

workflow {
    def paths = resolve_input_paths()

    def primary_run_name = params.run_name ?: (paths.bcl_dir ? run_name_for(file(paths.bcl_dir)) : 'run')

    def run_from = resolve_run_from()
    def extras   = resolve_extras()
    preflight(paths, run_from, extras)

    def si_indexes_fallback = load_si_indexes(projectDir.toString(), params.sequencer)
    if (run_from != 'bcl')
        log.info "no BCL dir, using --sequencer '${params.sequencer}' for i5 orientation"

    def all_rows = ([paths.samplesheet] + paths.extra_samplesheets).collectMany { ss_path ->
        parse_samplesheet(ss_path, si_indexes_fallback, primary_run_name, paths.adt_files_dir)
    }
    def conflicts = find_meta_conflicts(all_rows)
    if (conflicts)
        error "samplesheets disagree about the same library:\n  " + conflicts.join("\n  ")
    // A library listed in several samplesheets is counted once, from the first sheet's row.
    ch_meta = channel.fromList(all_rows.unique(false) { m -> m.id })

    if (run_from == 'cellranger') {
        ch_gex_outs = ch_meta
            .filter { meta -> counts_with_multi(meta) }
            .map { meta -> [meta.library_id, meta] }
            .groupTuple()
            .map { lid, metas -> [lid, metas, file("${paths.outs_dir}/${lid}/outs")] }
            .filter { _lid, _metas, outs -> outs.exists() }
            .ifEmpty { error "no cellranger outs found under ${paths.outs_dir}, expected ${paths.outs_dir}/<library_id>/outs" }

        ch_atac_outs = ch_meta
            .filter { meta -> meta.modality == 'ATAC' }
            .map { meta -> [meta, file("${paths.outs_dir}/${meta.library_id}_ATAC/outs")] }
            .branch { _meta, outs ->
                found:   outs.exists()
                missing: true
            }

        ch_atac_outs.missing.subscribe { meta, outs ->
            log.warn "no cellranger-atac outs for '${meta.library_id}' at ${outs}, skipping ATAC QC"
        }

        QC_GEX(ch_gex_outs)
        QC_ATAC(ch_atac_outs.found)

    } else {
        if (run_from == 'fastq') {
            ch_fastqs = channel.fromList(paths.fastq_dir)
                .combine(ch_meta)
                .map { fastq_dir, meta ->
                    def fqs = (new File(fastq_dir).listFiles() ?: []).findAll { f -> f.isFile() && f.name.endsWith('.fastq.gz') }
                    [meta, fastq_dir, sample_fastqs(meta.id, fqs).collect { f -> f.toPath() }.toSorted { p -> p.name }]
                }
                .filter { _meta, _fastq_dir, fqs -> !fqs.isEmpty() }
                .ifEmpty { error "no FASTQs matched any library under ${paths.fastq_dir}, expected files named <sample_id>_S<n>_*.fastq.gz" }
        } else {
            // Without --extra_samplesheets every extra flowcell is demultiplexed with the main samplesheet.
            def bcl_dirs   = [paths.bcl_dir] + paths.extra_bcl_dirs
            def bcl_sheets = [paths.samplesheet] + (paths.extra_samplesheets ?: [paths.samplesheet] * paths.extra_bcl_dirs.size())

            def meta_bcl_pairs = [bcl_dirs, bcl_sheets].transpose().collectMany { bcl_path, ss_path ->
                def flowcell_dir = file(bcl_path)
                def bcl_si       = load_si_indexes(projectDir.toString(), detect_sequencer(bcl_path, params.sequencer))
                parse_samplesheet(ss_path, bcl_si, run_name_for(flowcell_dir), paths.adt_files_dir).collect { meta -> [meta, flowcell_dir] }
            }

            // Each flowcell names its own FASTQs; counting output goes under the primary run name.
            ch_fastqs = DEMUX(channel.fromList(meta_bcl_pairs))
                .map { meta, fastq_dir, fqs -> [meta + [run_name: primary_run_name], fastq_dir, fqs] }
        }

        ch_routed = ch_fastqs
            .branch { meta, _fastq_dir, _fqs ->
                gex:      counts_with_multi(meta)
                atac:     meta.modality == 'ATAC'
                asap_adt: meta.assay == 'ASAP' && meta.modality in ['ADT', 'HTO']
                skip:     true
            }

        ch_routed.skip.subscribe { meta, _fastq_dir, _fqs ->
            log.warn "'${meta.id}' (assay=${meta.assay}, modality=${meta.modality}) matched no counting route, not counted"
        }

        ch_gex_libraries = by_library(ch_routed.gex)
            .map { lid, entries ->
                def by_mod  = entries.groupBy { e -> e.meta.modality }
                def metas   = by_mod.values().collect { es -> es[0].meta }.toSorted { m -> m.id }
                def adt_csv = metas.collect { m -> m.adt_csv_path }.find { p -> p }
                def fqs     = ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'].collect { mod ->
                    by_mod[mod]?.collectMany { e -> e.files } ?: [file('NO_FILE')]
                }
                [lid, metas, file(adt_csv ?: 'NO_FILE'), fqs]
            }

        ch_atac_libraries = by_library(ch_routed.atac)
            .map { _lid, entries -> [entries[0].meta, entries.collectMany { e -> e.files }] }

        ch_asap_adt = by_library(ch_routed.asap_adt)
            .map { lid, entries -> [lid, entries[0].meta, entries.collectMany { e -> e.files }.toSorted { f -> f.name }] }

        COUNT_GEX(ch_gex_libraries)
        COUNT_ATAC(ch_atac_libraries)

        // ADT is counted only for ASAP libraries whose ATAC half counted.
        COUNT_ADT(
            COUNT_ATAC.out
                .filter { meta, _outs -> meta.assay == 'ASAP' }
                .map { meta, _outs -> [meta.library_id] }
                .join(ch_asap_adt)
                .map { _lid, meta, fqs -> [meta, fqs] }
        )

        QC_GEX(COUNT_GEX.out)
        QC_ATAC(COUNT_ATAC.out)

        QUANT_EXTRA(COUNT_GEX.out, ch_routed.gex, QC_GEX.out, extras)
    }
}
