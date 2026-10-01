def load_si_indexes(String projectDir, String sequencer) {
    def i5_col   = sequencer == 'novaseq_x' ? 2 : 3
    def kit_rows = { String kit ->
        new File("${projectDir}/assets/indexes/${kit}_Set_A.csv").readLines()
            .findAll { line -> line.trim() && !line.startsWith('index_name') }
            .collect { line -> line.trim().split(',') }
    }
    [
        single: ['Single_Index_Kit_GA', 'Single_Index_Kit_NA']
            .collectMany(kit_rows).collectEntries { p -> [p[0], p[1..-1]] },
        dual:   ['Dual_Index_Kit_TT', 'Dual_Index_Kit_TN', 'Dual_Index_Kit_TS']
            .collectMany(kit_rows).collectEntries { p -> [p[0], [i7: p[1], i5: p[i5_col]]] },
    ]
}

def resolve_index(String index, Map si_indexes) {
    if (!index || index.toUpperCase() == 'NA')
        return [is_dual: false, rows: []]
    if (si_indexes.single.containsKey(index))
        return [is_dual: false, rows: si_indexes.single[index].collect { seq -> [i7: seq] }]
    if (si_indexes.dual.containsKey(index))
        return [is_dual: true, rows: [si_indexes.dual[index]]]
    [is_dual: false, rows: [[i7: index]]]
}

def detect_sequencer(String bcl_path, String fallback) {
    def runinfo    = new File("${bcl_path}/RunInfo.xml")
    def m          = runinfo.exists() ? runinfo.text =~ /<Instrument>([^<]+)<\/Instrument>/ : null
    def instrument = m ? m[0][1].trim() : ''
    def sequencer  = ['VH', 'NDX', 'LH', 'FS'].any { p -> instrument.startsWith(p) } ? 'novaseq_x'
                   : ['A', 'NB', 'NS', 'MN'].any { p -> instrument.startsWith(p) }  ? 'novaseq6000'
                   : null
    if (!sequencer) {
        log.warn "no known instrument ID in ${runinfo} ('${instrument}'), using --sequencer ${fallback}"
        return fallback
    }
    log.info "flow cell '${instrument}' is a ${sequencer} (i5 ${sequencer == 'novaseq_x' ? 'forward' : 'reverse-complement'})"
    sequencer
}

// <run>_bcl and <run>_fastq both name the run <run>.
def run_name_for(dir) {
    dir.name.replaceAll(/_(bcl|fastq)$/, '')
}

def fastq_dir_for(bcl_dir) {
    "${bcl_dir.parent}/${run_name_for(bcl_dir)}_fastq".toString()
}

// bcl-convert names files <Sample_ID>_S<n>_..., so a bare prefix would let lib1 claim lib10's FASTQs.
def sample_fastqs(String sample_id, files) {
    def pattern = java.util.regex.Pattern.compile("^" + java.util.regex.Pattern.quote(sample_id) + "_S\\d+_")
    files.findAll { f -> pattern.matcher(f.name).find() }
}
