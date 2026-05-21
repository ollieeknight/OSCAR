#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Imports ──────────────────────────────────────────────────────────────────
include { MULTIQC }        from './modules/demux'
include { DEMUX }          from './subworkflows/demux'
include { COUNT_GEX }      from './subworkflows/count_gex'
include { COUNT_ATAC }     from './subworkflows/count_atac'
include { COUNT_ADT }      from './subworkflows/count_adt'
include { QC_GEX }         from './subworkflows/qc_gex'
include { QC_ATAC }        from './subworkflows/qc_atac'

// ─── Helper functions ─────────────────────────────────────────────────────────

// Load 10x Genomics SI/DI index kit sequences from assets/indexes/ CSVs.
// Returns [single: {code: [seq1..4]}, dual: {code: {i7, i5}}]
// i5 column: NS6000 uses col3 (RC), NovaSeq X uses col4 (fwd)
def load_si_indexes(String projectDir, String sequencer) {
    def assets = "${projectDir}/assets/indexes"
    def si     = [single: [:], dual: [:]]
    def i5_col = (sequencer == 'novaseq_x') ? 2 : 3   // 1-indexed after name col

    ['Single_Index_Kit_T_Set_A.csv', 'Single_Index_Kit_N_Set_A.csv'].each { fname ->
        def f = new File("${assets}/${fname}")
        if (f.exists()) f.eachLine { line ->
            if (line.trim()) {
                def p = line.trim().split(',')
                si.single[p[0]] = p[1..-1]   // [seq1, seq2, seq3, seq4]
            }
        }
    }

    ['Dual_Index_Kit_TT_Set_A.csv', 'Dual_Index_Kit_Set_A.csv'].each { fname ->
        def f = new File("${assets}/${fname}")
        if (f.exists()) f.eachLine { line ->
            if (line.trim()) {
                def p = line.trim().split(',')
                si.dual[p[0]] = [i7: p[1], i5: p[i5_col]]
            }
        }
    }
    si
}

// Resolve a raw index value to BCL Convert rows.
// Returns [is_dual: bool, rows: [[i7: seq] or [i7: seq, i5: seq]]]
def resolve_index(String index, Map si_indexes) {
    if (si_indexes.single.containsKey(index))
        return [is_dual: false, rows: si_indexes.single[index].collect { [i7: it] }]
    if (si_indexes.dual.containsKey(index)) {
        def d = si_indexes.dual[index]
        return [is_dual: true, rows: [[i7: d.i7, i5: d.i5]]]
    }
    // Direct sequence (raw 8-mer or longer): single-index, i7 only
    return [is_dual: false, rows: [[i7: index]]]
}

def preflight_check() {
    if (!params.samplesheet) error "ERROR: --samplesheet is required"
    if (!file(params.samplesheet).exists()) error "ERROR: samplesheet not found: ${params.samplesheet}"

    if (params.from_fastq && params.from_cellranger)
        error "ERROR: --from-fastq and --from-cellranger are mutually exclusive"

    if (!params.from_fastq && !params.from_cellranger) {
        if (!params.bcl_dir) error "ERROR: --bcl_dir is required (or use --from-fastq / --from-cellranger)"
        if (!file(params.bcl_dir).exists()) error "ERROR: bcl_dir not found: ${params.bcl_dir}"
    }
    if (params.from_fastq && !params.fastq_dir)
        error "ERROR: --fastq_dir is required when --from-fastq is set"
    if (params.from_cellranger && !params.outs_dir)
        error "ERROR: --outs_dir is required when --from-cellranger is set"

    if (params.run_until && !['FASTQ', 'cellranger'].contains(params.run_until))
        error "ERROR: --run-until must be 'FASTQ' or 'cellranger' (got: '${params.run_until}')"
    if (params.from_fastq && params.run_until == 'FASTQ')
        error "ERROR: --run-until FASTQ has no effect when starting from FASTQs (--from-fastq)"
    if (params.from_cellranger && params.run_until)
        log.warn "WARNING: --run-until ignored — --from-cellranger starts at QC"
}

def preflight_samplesheet(String path) {
    def required = ['assay', 'experiment_id', 'historical_number', 'replicate',
                    'modality', 'chemistry', 'index_type', 'index',
                    'species', 'n_donors', 'adt_file']
    def valid_assays     = ['GEX', 'CITE', 'DOGMA', 'ATAC', 'Multiome', 'ASAP']
    def valid_modalities = ['GEX', 'ATAC', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR', 'GENO']

    def lines = new File(path).readLines()
    if (lines.isEmpty()) error "ERROR: samplesheet is empty: ${path}"

    def headers = lines[0].split(',').collect { it.trim() }
    def missing = required - headers
    if (missing) error "ERROR: samplesheet missing columns: ${missing.join(', ')}"

    lines.tail().eachWithIndex { line, i ->
        if (line.trim().isEmpty()) return
        def vals = line.split(',', -1).collect { it.trim() }
        if (vals.size() != headers.size())
            error "ERROR: samplesheet row ${i + 2}: expected ${headers.size()} fields, got ${vals.size()}"
        def row = [headers, vals].transpose().collectEntries()

        if (!valid_assays.any { it.equalsIgnoreCase(row.assay) })
            error "ERROR: row ${i + 2}: unknown assay '${row.assay}'. Valid: ${valid_assays.join(', ')}"
        if (!valid_modalities.contains(row.modality))
            error "ERROR: row ${i + 2}: unknown modality '${row.modality}'. Valid: ${valid_modalities.join(', ')}"
        if (!['human', 'mouse'].any { it.equalsIgnoreCase(row.species) })
            error "ERROR: row ${i + 2}: unknown species '${row.species}'. Valid: human, mouse"
    }
}

def parse_row(row, Map si_indexes) {
    def n_donors = (row.n_donors == null || row.n_donors.trim() in ['NA', '', 'na']) \
        ? 1 : row.n_donors.trim().toInteger()
    def index    = row.index.trim()
    def adt_file = row.adt_file?.trim() ?: null
    def adt_csv_path = (adt_file && params.adt_files_dir) \
        ? file("${params.adt_files_dir}/${adt_file}.csv").toAbsolutePath().toString() : null
    [
        id:               "${row.assay}_${row.experiment_id}_exp${row.historical_number}_lib${row.replicate}_${row.modality}",
        library_id:       "${row.assay}_${row.experiment_id}_exp${row.historical_number}_lib${row.replicate}",
        assay:            row.assay.trim(),
        experiment_id:    row.experiment_id.trim(),
        historical_number: row.historical_number.trim(),
        replicate:        row.replicate.trim(),
        modality:         row.modality.trim(),
        chemistry:        row.chemistry.trim(),
        index_type:       row.index_type.trim(),
        index:            index,
        index_seqs:       resolve_index(index, si_indexes),
        species:          row.species.trim().toLowerCase(),
        n_donors:         n_donors,
        adt_file:         adt_file,
        adt_csv_path:     adt_csv_path
    ]
}

// ─── Sequencer auto-detection ────────────────────────────────────────────────
// Reads <Instrument> from RunInfo.xml and maps the prefix to the i5 orientation
// used by load_si_indexes().
//
// Instrument ID prefixes:
//   VH  → NovaSeq X / X Plus   → i5 forward  → 'novaseq_x'
//   A   → NovaSeq 6000          → i5 RC        → 'novaseq6000'
//   LH  → NovaSeq 6000          → i5 RC        → 'novaseq6000'
//   MN  → MiniSeq               → i5 RC        → 'novaseq6000' (fallback)
//   NB  → NextSeq 550           → i5 RC        → 'novaseq6000'
//   NS  → NextSeq 500           → i5 RC        → 'novaseq6000'
//   NDX → NextSeq 2000/1000     → i5 forward   → 'novaseq_x'
// If RunInfo.xml is absent or unparseable, warns and uses params.sequencer.

def detect_sequencer(String bcl_path) {
    def runinfo = new File("${bcl_path}/RunInfo.xml")
    if (!runinfo.exists()) {
        log.warn "WARNING: RunInfo.xml not found in ${bcl_path}; falling back to params.sequencer='${params.sequencer}'"
        return params.sequencer
    }
    def text = runinfo.text
    def m    = text =~ /<Instrument>([^<]+)<\/Instrument>/
    if (!m) {
        log.warn "WARNING: <Instrument> tag not found in ${bcl_path}/RunInfo.xml; falling back to params.sequencer='${params.sequencer}'"
        return params.sequencer
    }
    def instrument_id = m[0][1].trim()
    def sequencer
    if      (instrument_id.startsWith('VH'))  sequencer = 'novaseq_x'   // NovaSeq X / X Plus
    else if (instrument_id.startsWith('NDX')) sequencer = 'novaseq_x'   // NextSeq 2000/1000
    else if (instrument_id.startsWith('A'))   sequencer = 'novaseq6000' // NovaSeq 6000
    else if (instrument_id.startsWith('LH'))  sequencer = 'novaseq6000' // NovaSeq 6000 (XLEAP)
    else if (instrument_id.startsWith('NB'))  sequencer = 'novaseq6000' // NextSeq 550
    else if (instrument_id.startsWith('NS'))  sequencer = 'novaseq6000' // NextSeq 500
    else if (instrument_id.startsWith('MN'))  sequencer = 'novaseq6000' // MiniSeq
    else {
        log.warn "WARNING: Unrecognised instrument ID '${instrument_id}' in ${bcl_path}/RunInfo.xml; falling back to params.sequencer='${params.sequencer}'"
        sequencer = params.sequencer
    }
    log.info "INFO: Detected instrument '${instrument_id}' → sequencer mode '${sequencer}' (i5 ${sequencer == 'novaseq_x' ? 'forward' : 'reverse-complement'})"
    return sequencer
}

// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    preflight_check()

    def all_ss_paths = [params.samplesheet]
    if (params.extra_samplesheets)
        all_ss_paths += params.extra_samplesheets.split(',').collect { it.trim() }
    all_ss_paths.each { preflight_samplesheet(it) }

    // ── Resolve sequencer / i5 orientation (non-BCL fallback) ────────────────
    // In BCL mode, sequencer is auto-detected per BCL dir inside the BCL branch
    // below (each flowcell may come from a different instrument). For --from-fastq
    // and --from-cellranger there is no BCL dir, so we fall back to params.sequencer.
    def si_indexes_fallback = load_si_indexes(projectDir.toString(), params.sequencer)
    if (params.from_fastq || params.from_cellranger)
        log.info "INFO: No BCL dir available — using params.sequencer='${params.sequencer}' for i5 orientation"

    // ── Parse samplesheets → ch_meta ─────────────────────────────────────────
    // For from_fastq / from_cellranger, ch_meta is set here and used downstream.
    // For BCL mode, the BCL branch below rebuilds ch_meta with per-instrument
    // index sequences — this initial set is overridden there.
    def _all_rows = []
    all_ss_paths.each { ss_path ->
        def lines = new File(ss_path).readLines()
        if (!lines.isEmpty()) {
            def hdrs = lines[0].split(',').collect { it.trim() }
            lines.tail().each { line ->
                if (!line.trim().isEmpty()) {
                    def vals = line.split(',', -1).collect { it.trim() }
                    _all_rows << parse_row([hdrs, vals].transpose().collectEntries(), si_indexes_fallback)
                }
            }
        }
    }
    Channel.fromList(_all_rows).set { ch_meta }

    // ── Entry point: --from-cellranger (QC only, run_until ignored) ───────────

    if (params.from_cellranger) {
        ch_meta
            .filter { meta ->
                meta.modality in ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] \
                    && meta.assay != 'ASAP'
            }
            .map { meta -> [meta.library_id, meta] }
            .groupTuple(by: 0)
            .map { lid, metas -> [lid, metas, file("${params.outs_dir}/${lid}/outs")] }
            .filter { lid, metas, outs -> outs.exists() }
            .set { ch_gex_outs }

        ch_meta
            .filter { meta -> meta.modality == 'ATAC' }
            .map { meta -> [meta, file("${params.outs_dir}/${meta.library_id}_ATAC/outs")] }
            .filter { meta, outs -> outs.exists() }
            .set { ch_atac_outs }

        QC_GEX(ch_gex_outs)
        QC_ATAC(ch_atac_outs)

        QC_GEX.out.logs.flatten().collect().set { ch_multiqc_input }
        MULTIQC(ch_multiqc_input)

    } else {
        // ── Obtain FASTQs (BCL demux or pre-existing) ─────────────────────────
        if (params.from_fastq) {
            ch_meta
                .map { meta ->
                    def fqs = file("${params.fastq_dir}/**/${meta.id}*.fastq.gz")
                    fqs = fqs instanceof List ? fqs : (fqs.exists() ? [fqs] : [])
                    def parents = fqs.collect { it.parent.toString() }.unique()
                    parents.collect { pdir ->
                        def matched_fqs = fqs.findAll { it.parent.toString() == pdir }
                        [meta, pdir, matched_fqs]
                    }
                }
                .flatMap()
                .filter { meta, fastq_dir, fqs -> !fqs.isEmpty() }
                .set { ch_fastqs }
        } else {
            def bcl_paths = [params.bcl_dir]
            def bcl_ss    = [params.samplesheet]
            if (params.extra_bcl_dirs) {
                def extra_bcls = params.extra_bcl_dirs.split(',').collect { it.trim() }
                bcl_paths += extra_bcls
                if (params.extra_samplesheets) {
                    def extra_sss = params.extra_samplesheets.split(',').collect { it.trim() }
                    if (extra_sss.size() != extra_bcls.size())
                        error "ERROR: --extra_samplesheets count (${extra_sss.size()}) must match --extra_bcl_dirs (${extra_bcls.size()})"
                    bcl_ss += extra_sss
                } else {
                    bcl_ss += extra_bcls.collect { params.samplesheet }
                }
            }

            def _meta_bcl_pairs = []
            def _bcl_rows       = []   // used to rebuild ch_meta with correctly-resolved index_seqs
            [bcl_paths, bcl_ss].transpose().each { bcl_path, ss_path ->
                def bcl_dir       = file(bcl_path)
                // Detect sequencer per BCL dir — each flowcell may originate from a
                // different instrument (e.g. mixing NovaSeq X and NovaSeq 6000 runs).
                def bcl_si        = load_si_indexes(projectDir.toString(), detect_sequencer(bcl_path))
                def lines         = new File(ss_path).readLines()
                if (!lines.isEmpty()) {
                    def hdrs = lines[0].split(',').collect { it.trim() }
                    lines.tail().each { line ->
                        if (!line.trim().isEmpty()) {
                            def vals = line.split(',', -1).collect { it.trim() }
                            def meta = parse_row([hdrs, vals].transpose().collectEntries(), bcl_si)
                            _meta_bcl_pairs << [meta, bcl_dir]
                            _bcl_rows       << meta
                        }
                    }
                }
            }
            // Override ch_meta with the correctly-detected index sequences from each
            // BCL dir. Deduplicate by meta.id (same sample listed in multiple flowcell
            // samplesheets) while preserving insertion order.
            def seen_ids    = [] as Set
            def unique_rows = _bcl_rows.findAll { seen_ids.add(it.id) }
            Channel.fromList(unique_rows).set { ch_meta }
            Channel.fromList(_meta_bcl_pairs).set { ch_meta_bcl }

            DEMUX(ch_meta_bcl)
            ch_fastqs = DEMUX.out.fastqs
        }

        // ── --run-until FASTQ: stop after demux ───────────────────────────────
        if (params.run_until == 'FASTQ') {
            // Only BCL path produces falco reports; from_fastq has no demux step
            DEMUX.out.falco_reports.flatten().collect().set { ch_multiqc_input }
            MULTIQC(ch_multiqc_input)

        } else {
            // ── Count ─────────────────────────────────────────────────────────
            ch_fastqs
                .branch { meta, fastq_dir, fqs ->
                    gex:      meta.modality in ['GEX', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] \
                              && meta.assay != 'ASAP'
                    atac:     meta.modality == 'ATAC'
                    asap_adt: meta.assay == 'ASAP' && meta.modality in ['ADT', 'HTO']
                    skip:     true
                }
                .set { ch_routed }

            // GEX: group by library_id, deduplicate metas by modality, collect unique fastq dirs.
            // Each fastq dir is a published path string — no staging, so no filename collision
            // when the same library was sequenced on multiple flowcells.
            ch_routed.gex
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, fastq_dir] }
                .groupTuple(by: 0)
                .map { lid, metas, fastq_dirs ->
                    def seen         = [] as Set
                    def unique_metas = metas.findAll { m -> seen.add(m.modality) }
                    def unique_dirs  = fastq_dirs.toSet().toList()
                    def adt_csv_path = metas.collect { it.adt_csv_path }.find { it }
                    def adt_csv      = adt_csv_path ? file(adt_csv_path) : file('NO_FILE')
                    [lid, unique_metas, unique_dirs, adt_csv]
                }
                .set { ch_gex_libraries }

            // ATAC: group by library_id to collect dirs from multiple flowcells
            ch_routed.atac
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, fastq_dir] }
                .groupTuple(by: 0)
                .map { lid, metas, fastq_dirs -> [metas[0], fastq_dirs.toSet().toList()] }
                .set { ch_atac_libraries }

            COUNT_GEX(ch_gex_libraries)
            COUNT_ATAC(ch_atac_libraries)

            ch_asap_atac_outs = COUNT_ATAC.out.outs
                .filter { meta, outs -> meta.assay == 'ASAP' }
                .map    { meta, outs -> [meta.library_id, meta, outs] }

            ch_asap_adt_fastqs = ch_routed.asap_adt
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, fqs] }
                .groupTuple(by: 0)
                .map { lid, metas, fq_lists -> [lid, metas[0], fq_lists.flatten()] }

            COUNT_ADT(
                ch_asap_atac_outs
                    .join(ch_asap_adt_fastqs, by: 0, failOnDuplicate: false, failOnMismatch: false)
                    .map { lid, atac_meta, outs, adt_meta, adt_fqs ->
                        [atac_meta, outs, adt_meta, adt_fqs]
                    }
            )

            // ── --run-until cellranger: stop after counting ───────────────────
            if (params.run_until == 'cellranger') {
                if (!params.from_fastq) {
                    DEMUX.out.falco_reports.flatten().collect().set { ch_multiqc_input }
                    MULTIQC(ch_multiqc_input)
                }

            } else {
                // ── Full run: QC + MultiQC ────────────────────────────────────
                QC_GEX(COUNT_GEX.out.outs)
                QC_ATAC(COUNT_ATAC.out.outs)

                Channel.empty()
                    .mix(params.from_fastq ? Channel.empty() : DEMUX.out.falco_reports.flatten())
                    .mix(QC_GEX.out.logs.flatten())
                    .collect()
                    .set { ch_multiqc_input }

                MULTIQC(ch_multiqc_input)
            }
        }
    }
}
