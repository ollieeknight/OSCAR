def load_si_indexes(String projectDir, String sequencer) {
    def assets = "${projectDir}/assets/indexes"
    def si     = [single: [:], dual: [:]]
    def i5_col = (sequencer == 'novaseq_x') ? 2 : 3

    ['Single_Index_Kit_GA_Set_A.csv', 'Single_Index_Kit_NA_Set_A.csv'].each { fname ->
        def f = new File("${assets}/${fname}")
        if (f.exists()) f.eachLine { line ->
            if (line.trim() && !line.startsWith('index_name')) {
                def p = line.trim().split(',')
                si.single[p[0]] = p[1..-1]
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
    return [is_dual: false, rows: [[i7: index]]]
}

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
    if      (instrument_id.startsWith('VH'))  sequencer = 'novaseq_x'
    else if (instrument_id.startsWith('NDX')) sequencer = 'novaseq_x'
    else if (instrument_id.startsWith('A'))   sequencer = 'novaseq6000'
    else if (instrument_id.startsWith('LH'))  sequencer = 'novaseq_x'
    else if (instrument_id.startsWith('NB'))  sequencer = 'novaseq6000'
    else if (instrument_id.startsWith('NS'))  sequencer = 'novaseq6000'
    else if (instrument_id.startsWith('MN'))  sequencer = 'novaseq6000'
    else if (instrument_id.startsWith('FS'))  sequencer = 'novaseq_x'
    else {
        log.warn "WARNING: Unrecognised instrument ID '${instrument_id}' in ${bcl_path}/RunInfo.xml; falling back to params.sequencer='${fallback}'"
        sequencer = fallback
    }
    log.info "INFO: Flow cell '${instrument_id}' indicates a '${sequencer}' (i5 set to ${sequencer == 'novaseq_x' ? 'forward' : 'reverse-complement'})"
    return sequencer
}
