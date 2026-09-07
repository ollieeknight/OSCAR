// ─── Index kits and sequencer detection ──────────────────────────────────────
// Loading of 10x SI/DI index kit sequences, resolution of a samplesheet index
// value to BCL Convert rows, and instrument-ID → i5-orientation detection.

// Load 10x Genomics SI/DI index kit sequences from assets/indexes/ CSVs.
// Returns [single: {code: [seq1..4]}, dual: {code: {i7, i5}}]
//
// i5 column: 10x ships two i5 orientations per dual-index kit.
//   col 3 (index 2) = workflow A, forward   → NovaSeq X, NextSeq 1000/2000, MiSeq
//   col 4 (index 3) = workflow B, RC        → NovaSeq 6000, NextSeq 500/550
//
// Some kit CSVs (TS) carry a header row and some (TT, TN, GA, NA) do not, so
// header lines are skipped explicitly rather than by position.
def load_si_indexes(String projectDir, String sequencer) {
    def assets = "${projectDir}/assets/indexes"
    def si     = [single: [:], dual: [:]]
    def i5_col = (sequencer == 'novaseq_x') ? 2 : 3   // 1-indexed after name col

    ['Single_Index_Kit_GA_Set_A.csv', 'Single_Index_Kit_NA_Set_A.csv'].each { fname ->
        def f = new File("${assets}/${fname}")
        if (f.exists()) f.eachLine { line ->
            if (line.trim() && !line.startsWith('index_name')) {
                def p = line.trim().split(',')
                si.single[p[0]] = p[1..-1]   // [seq1, seq2, seq3, seq4]
            }
        }
    }

    ['Dual_Index_Kit_TT_Set_A.csv', 'Dual_Index_Kit_TN_Set_A.csv', 'Dual_Index_Kit_TS_Set_A.csv'].each { fname ->
        def f = new File("${assets}/${fname}")
        if (f.exists()) f.eachLine { line ->
            if (line.trim() && !line.startsWith('index_name')) {
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
    if (index == null || index.toUpperCase() == 'NA' || index == '') {
        return [is_dual: false, rows: []]
    }
    if (si_indexes.single.containsKey(index))
        return [is_dual: false, rows: si_indexes.single[index].collect { seq -> [i7: seq] }]
    if (si_indexes.dual.containsKey(index)) {
        def d = si_indexes.dual[index]
        return [is_dual: true, rows: [[i7: d.i7, i5: d.i5]]]
    }
    // Direct sequence (raw 8-mer or longer): single-index, i7 only
    return [is_dual: false, rows: [[i7: index]]]
}

// ─── Sequencer auto-detection ────────────────────────────────────────────────
// Reads <Instrument> from RunInfo.xml and maps the prefix to the i5 orientation
// used by load_si_indexes().
//
// Instrument ID prefixes:
//   VH  → NovaSeq X / X Plus   → i5 forward  → 'novaseq_x'
//   A   → NovaSeq 6000          → i5 RC        → 'novaseq6000'
//   LH  → NovaSeq X LEAP         → i5 forward   → 'novaseq_x'
//   MN  → MiniSeq               → i5 RC        → 'novaseq6000' (fallback)
//   NB  → NextSeq 550           → i5 RC        → 'novaseq6000'
//   NS  → NextSeq 500           → i5 RC        → 'novaseq6000'
//   NDX → NextSeq 2000/1000     → i5 forward   → 'novaseq_x'
// If RunInfo.xml is absent or unparseable, warns and uses params.sequencer.
def detect_sequencer(String bcl_path, String fallback) {
    def runinfo = new File("${bcl_path}/RunInfo.xml")
    if (!runinfo.exists()) {
        log.warn "WARNING: RunInfo.xml not found in ${bcl_path}; falling back to params.sequencer='${fallback}'"
        return fallback
    }
    def text = runinfo.text
    def m    = text =~ /<Instrument>([^<]+)<\/Instrument>/
    if (!m) {
        log.warn "WARNING: <Instrument> tag not found in ${bcl_path}/RunInfo.xml; falling back to params.sequencer='${fallback}'"
        return fallback
    }
    def instrument_id = m[0][1].trim()
    def sequencer
    if      (instrument_id.startsWith('VH'))  sequencer = 'novaseq_x'   // NovaSeq X / X Plus
    else if (instrument_id.startsWith('NDX')) sequencer = 'novaseq_x'   // NextSeq 2000/1000
    else if (instrument_id.startsWith('A'))   sequencer = 'novaseq6000' // NovaSeq 6000
    else if (instrument_id.startsWith('LH'))  sequencer = 'novaseq_x'   // NovaSeq X LEAP (i5 forward; BCL Convert auto-RCs via IsReverseComplement)
    else if (instrument_id.startsWith('NB'))  sequencer = 'novaseq6000' // NextSeq 550
    else if (instrument_id.startsWith('NS'))  sequencer = 'novaseq6000' // NextSeq 500
    else if (instrument_id.startsWith('MN'))  sequencer = 'novaseq6000' // MiniSeq
    else if (instrument_id.startsWith('FS'))  sequencer = 'novaseq_x'   // iSeq 100
    else {
        log.warn "WARNING: Unrecognised instrument ID '${instrument_id}' in ${bcl_path}/RunInfo.xml; falling back to params.sequencer='${fallback}'"
        sequencer = fallback
    }
    log.info "INFO: Detected instrument '${instrument_id}' → sequencer mode '${sequencer}' (i5 ${sequencer == 'novaseq_x' ? 'forward' : 'reverse-complement'})"
    return sequencer
}
