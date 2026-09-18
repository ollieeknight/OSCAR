#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEMUX }          from './subworkflows/demux'
include { COUNT_GEX }      from './subworkflows/count_gex'
include { COUNT_ATAC }     from './subworkflows/count_atac'
include { COUNT_ADT }      from './subworkflows/count_adt'
include { QC_GEX }         from './subworkflows/qc_gex'
include { QC_ATAC }        from './subworkflows/qc_atac'

include { QUANT_EXTRA }    from './subworkflows/quant_extra'

include { load_si_indexes; detect_sequencer } from './lib/indexes'
include { preflight_samplesheet; parse_samplesheet; find_meta_conflicts } from './lib/samplesheet'


def preflight_check(Map paths) {
    if (!paths.samplesheet) error "ERROR: --samplesheet is required"
    if (!file(paths.samplesheet).exists()) error "ERROR: samplesheet not found: ${paths.samplesheet}"

    def run_from = resolve_run_from()
    if (run_from == 'bcl') {
        if (!paths.bcl_dir) error "ERROR: --bcl_dir is required when --run_from bcl"
        if (!file(paths.bcl_dir).exists()) error "ERROR: bcl_dir not found: ${paths.bcl_dir}"
    }
    if (run_from == 'fastq') {
        if (!paths.fastq_dir) error "ERROR: --fastq_dir is required when --run_from fastq"
        paths.fastq_dir.each { d ->
            if (!file(d).exists()) error "ERROR: fastq_dir not found: ${d}"
        }
    }
    if (run_from == 'cellranger' && !paths.outs_dir)
        error "ERROR: --outs_dir is required when --run_from cellranger"

    def extras = resolve_extras()
    if ('viral' in extras) {
        if (!params.viral_piscem_index || !file(params.viral_piscem_index).exists())
            error "ERROR: --extras viral requires viral_piscem_index: ${params.viral_piscem_index}"
        if (!params.viral_t2g || !file(params.viral_t2g).exists())
            error "ERROR: --extras viral requires viral_t2g: ${params.viral_t2g}"
        if (!(params.bamtofastq_bin ==~ /^https?:.*/) && !file(params.bamtofastq_bin).exists())
            error "ERROR: --extras viral requires bamtofastq_bin: ${params.bamtofastq_bin}"
    }
    if ('velocity' in extras) {
        def spliceu = [human: params.spliceu_index_human, mouse: params.spliceu_index_mouse]
        def missing = spliceu.findAll { _sp, idx -> !idx || !file(idx).exists() }
        if (missing.size() == spliceu.size())
            error "ERROR: --extras velocity requires a spliceu index: ${spliceu.values().join(', ')}"
        missing.each { sp, idx -> log.warn "WARNING: --extras velocity: no ${sp} spliceu index (${idx}) — ${sp} libraries will fail" }
    }
}


def resolve_run_from() {
    def rf = (params.run_from ?: 'bcl').toString().toLowerCase()
    if (!['bcl', 'fastq', 'cellranger'].contains(rf))
        error "ERROR: --run_from must be 'bcl', 'fastq' or 'cellranger' (got: '${params.run_from}')"
    return rf
}

def resolve_extras() {
    def extras = (params.extras ?: '').toString().tokenize(',')*.trim()*.toLowerCase().findAll { e -> e }
    def unknown = extras - ['velocity', 'viral']
    if (unknown)
        error "ERROR: unknown --extras: ${unknown.join(', ')} (valid: velocity, viral)"
    return extras
}

def preflight_flex(String ss_path) {
    def lines = new File(ss_path).readLines().findAll { line -> !line.trim().isEmpty() }
    def headers = lines[0].split(',').collect { h -> h.trim() }
    def has_flex = lines.tail().any { line ->
        def vals = line.split(',', -1)
        def row = [headers, vals].transpose().collectEntries()
        row.assay?.trim()?.equalsIgnoreCase('Flex')
    }
    if (!has_flex) return

    if (!(params.flex_backend in ['cellranger', 'cyto', 'both']))
        error "ERROR: --flex_backend must be 'cellranger', 'cyto', or 'both' (got '${params.flex_backend}')"

    def uses_cr   = params.flex_backend in ['cellranger', 'both']
    def uses_cyto = params.flex_backend in ['cyto', 'both']

    if (uses_cr) {
        if (!params.flex_probe_set)
            error "ERROR: Flex run requires --flex_probe_set (standard 10x probe set CSV) when flex_backend includes cellranger."
        if (!file(params.flex_probe_set).exists())
            error "ERROR: --flex_probe_set not found: ${params.flex_probe_set}"
    }
    if (params.flex_probe_set_custom && !file(params.flex_probe_set_custom).exists())
        error "ERROR: --flex_probe_set_custom not found: ${params.flex_probe_set_custom}"
    if (params.flex_samples_file && !file(params.flex_samples_file).exists())
        error "ERROR: --flex_samples_file not found: ${params.flex_samples_file}"

    if (uses_cyto) {
        if (!params.container_cyto || !file(params.container_cyto).exists())
            error "ERROR: cyto container not found at '${params.container_cyto}'. " +
                  "Build it: apptainer build ${params.container_cyto} OSCAR/containers/cyto.def"
    }
}

def refs_for(meta) {
    def is_human = meta.species == 'human'
    [
        gex: is_human ? params.ref_human     : params.ref_mouse,
        vdj: is_human ? params.ref_vdj_human : params.ref_vdj_mouse
    ]
}

def to_abs_path(pth) {
    pth ? file(pth).toAbsolutePath().toString() : null
}

def to_abs_list(csv) {
    csv ? csv.split(',').collect { pth -> file(pth.trim()).toAbsolutePath().toString() } : []
}

def resolve_input_paths() {
    [
        samplesheet:        to_abs_path(params.samplesheet),
        bcl_dir:            to_abs_path(params.bcl_dir),
        fastq_dir:          to_abs_list(params.fastq_dir),
        outs_dir:           to_abs_path(params.outs_dir),
        adt_files_dir:      to_abs_path(params.adt_files_dir),
        extra_samplesheets: to_abs_list(params.extra_samplesheets),
        extra_bcl_dirs:     to_abs_list(params.extra_bcl_dirs),
    ]
}


workflow {
    def paths = resolve_input_paths()

    def primary_run_name = params.run_name
    if (!primary_run_name || primary_run_name == 'null') {
        primary_run_name = paths.bcl_dir
            ? file(paths.bcl_dir).name.replaceAll(/_bcl$/, '')
            : 'run'
    }

    preflight_check(paths)

    def run_from = resolve_run_from()
    def extras   = resolve_extras()

    def all_ss_paths = [paths.samplesheet] + paths.extra_samplesheets
    all_ss_paths.each { ss -> preflight_samplesheet(ss) }
    all_ss_paths.each { ss -> preflight_flex(ss) }

    def si_indexes_fallback = load_si_indexes(projectDir.toString(), params.sequencer)
    if (run_from != 'bcl')
        log.info "INFO: No BCL dir available — using params.sequencer='${params.sequencer}' for i5 orientation"

    def all_rows = all_ss_paths.collectMany { ss_path ->
        parse_samplesheet(ss_path, si_indexes_fallback, primary_run_name, paths.adt_files_dir)
    }
    channel.fromList(all_rows).set { ch_meta }


    if (run_from == 'cellranger') {
        ch_meta
            .filter { meta ->
                (meta.modality in ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] \
                    && meta.assay != 'ASAP') || meta.assay == 'Flex'
            }
            .map { meta -> [meta.library_id, meta] }
            .groupTuple(by: 0)
            .map { lid, metas -> [lid, metas, file("${paths.outs_dir}/${lid}/outs")] }
            .filter { _lid, _metas, outs -> outs.exists() }
            .ifEmpty { error "ERROR: no cellranger outs found under ${paths.outs_dir} — expected ${paths.outs_dir}/<library_id>/outs" }
            .set { ch_gex_outs }

        ch_meta
            .filter { meta -> meta.modality == 'ATAC' }
            .map { meta -> [meta, file("${paths.outs_dir}/${meta.library_id}_ATAC/outs")] }
            .branch { _meta, outs ->
                found:   outs.exists()
                missing: true
            }
            .set { ch_atac_routed }

        ch_atac_routed.missing.subscribe { meta, outs ->
            log.warn "WARN: no cellranger-atac outs for '${meta.library_id}' at ${outs} — skipping ATAC QC"
        }

        ch_atac_routed.found.set { ch_atac_outs }

        QC_GEX(ch_gex_outs)
        QC_ATAC(ch_atac_outs)

    } else {
        if (run_from == 'fastq') {
            channel.fromList(paths.fastq_dir)
                .combine(ch_meta)
                .map { fastq_dir, meta ->
                    def id_re = java.util.regex.Pattern.compile(
                        "^" + java.util.regex.Pattern.quote(meta.id) + "_S\\d+_")
                    def fqs = (new File(fastq_dir).listFiles() ?: [])
                        .findAll { f -> f.isFile() && f.name.endsWith('.fastq.gz') \
                                        && id_re.matcher(f.name).find() }
                        .collect { f -> f.toPath() }
                        .toSorted { p -> p.getFileName().toString() }
                    [meta, fastq_dir, fqs]
                }
                .filter { _meta, _fastq_dir, fqs -> !fqs.isEmpty() }
                .ifEmpty { error "ERROR: no FASTQs matched any library under ${paths.fastq_dir} — expected files named <sample_id>_S<n>_*.fastq.gz" }
                .set { ch_fastqs }
        } else {
            def bcl_paths = [paths.bcl_dir]
            def bcl_ss    = [paths.samplesheet]
            if (paths.extra_bcl_dirs) {
                bcl_paths += paths.extra_bcl_dirs
                if (paths.extra_samplesheets) {
                    if (paths.extra_samplesheets.size() != paths.extra_bcl_dirs.size())
                        error "ERROR: --extra_samplesheets count (${paths.extra_samplesheets.size()}) must match --extra_bcl_dirs (${paths.extra_bcl_dirs.size()})"
                    bcl_ss += paths.extra_samplesheets
                } else {
                    bcl_ss += paths.extra_bcl_dirs.collect { _b -> paths.samplesheet }
                }
            }

            def meta_bcl_pairs = []
            def bcl_rows       = []
            [bcl_paths, bcl_ss].transpose().each { bcl_path, ss_path ->
                def flowcell_dir = file(bcl_path)
                def bcl_si  = load_si_indexes(projectDir.toString(), detect_sequencer(bcl_path, params.sequencer))
                def run_nm  = flowcell_dir.name.replaceAll(/_bcl$/, '')
                parse_samplesheet(ss_path, bcl_si, run_nm, paths.adt_files_dir).each { meta ->
                    meta_bcl_pairs << [meta, flowcell_dir]
                    bcl_rows       << meta
                }
            }
            def conflicts = find_meta_conflicts(bcl_rows)
            if (conflicts)
                error "ERROR: samplesheets disagree about the same library:\n  " + conflicts.join("\n  ")

            def seen_ids    = [] as Set
            def unique_rows = bcl_rows.findAll { m -> seen_ids.add(m.id) }
            channel.fromList(unique_rows).set { ch_meta }
            channel.fromList(meta_bcl_pairs).set { ch_meta_bcl }

            DEMUX(ch_meta_bcl)
            ch_fastqs = DEMUX.out.fastqs
        }

        ch_fastqs
            .branch { meta, _fastq_dir, _fqs ->
                gex:      (meta.modality in ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] \
                          && meta.assay != 'ASAP') || meta.assay == 'Flex'
                atac:     meta.modality == 'ATAC'
                asap_adt: meta.assay == 'ASAP' && meta.modality in ['ADT', 'HTO']
                skip:     true
            }
            .set { ch_routed }

        ch_routed.skip.subscribe { meta, _fastq_dir, _fqs ->
            log.warn "WARN: '${meta.id}' (assay=${meta.assay}, modality=${meta.modality}) matched no counting route — not counted"
        }

        ch_routed.gex
            .tap { ch_gex_for_velocity }
            .map { meta, _fastq_dir, fqs ->
                def files = (fqs instanceof List ? fqs : [fqs]).toSorted { f -> f.name }
                [meta.library_id, [modality: meta.modality, meta: meta, files: files]]
            }
            .groupTuple(by: 0)
            .map { lid, entries ->
                def mod_data = [:].withDefault { [meta: null, entries: [], files: []] }
                entries.each { e ->
                    if (!mod_data[e.modality].meta) mod_data[e.modality].meta = e.meta
                    mod_data[e.modality].entries << e
                }
                mod_data.each { _mod, d ->
                    d.entries = d.entries.toSorted { e -> e.files[0].toUriString() }
                    d.files   = d.entries.collectMany { e -> e.files }
                }

                def all_metas = []
                mod_data.each { _mod, d -> all_metas << d.meta }
                all_metas = all_metas.toSorted { m -> m.id }
                all_metas = all_metas.collect { m -> m + [run_name: primary_run_name] }

                def meta         = all_metas.find { m -> m.modality == 'GEX' } ?: all_metas[0]
                def adt_csv_path = all_metas.collect { m -> m.adt_csv_path }.find { p -> p }
                def adt_csv      = adt_csv_path ? file(adt_csv_path) : file('NO_FILE')

                [lid, all_metas, refs_for(meta), adt_csv,
                 (mod_data['GEX']?.files)    ?: [file('NO_FILE')],
                 (mod_data['ADT']?.files)    ?: [file('NO_FILE')],
                 (mod_data['HTO']?.files)    ?: [file('NO_FILE')],
                 (mod_data['VDJ-T']?.files)  ?: [file('NO_FILE')],
                 (mod_data['VDJ-B']?.files)  ?: [file('NO_FILE')],
                 (mod_data['CRISPR']?.files) ?: [file('NO_FILE')]]
            }
            .set { ch_gex_libraries }

        ch_routed.atac
            .map { meta, _fastq_dir, fqs ->
                def files = (fqs instanceof List ? fqs : [fqs]).toSorted { f -> f.name }
                [meta.library_id, [meta: meta, files: files]]
            }
            .groupTuple(by: 0)
            .map { _lid, entries ->
                def meta = entries[0].meta + [run_name: primary_run_name]

                def all_files = entries.toSorted { e -> e.files[0].toUriString() }.collectMany { e -> e.files }
                [meta, all_files]
            }
            .set { ch_atac_libraries }

        COUNT_GEX(ch_gex_libraries)
        COUNT_ATAC(ch_atac_libraries)

        ch_asap_atac_outs = COUNT_ATAC.out
            .filter { meta, _outs -> meta.assay == 'ASAP' }
            .map    { meta, outs -> [meta.library_id, meta, outs] }

        ch_asap_adt_fastqs = ch_routed.asap_adt
            .map { meta, _fastq_dirs, fqs -> [meta.library_id, meta, fqs] }
            .groupTuple(by: 0)
            .map { lid, metas, fq_lists -> [lid, metas[0], fq_lists.flatten()] }

        COUNT_ADT(
            ch_asap_atac_outs
                .join(ch_asap_adt_fastqs, by: 0, failOnDuplicate: false, failOnMismatch: false)
                .map { _lid, atac_meta, outs, adt_meta, adt_fqs ->
                    [atac_meta, outs, adt_meta, adt_fqs]
                }
        )

        QC_GEX(COUNT_GEX.out)
        QC_ATAC(COUNT_ATAC.out)

        QUANT_EXTRA(
            COUNT_GEX.out,
            ch_gex_for_velocity,
            QC_GEX.out.barcodes,
            extras,
            primary_run_name
        )
    }
}
