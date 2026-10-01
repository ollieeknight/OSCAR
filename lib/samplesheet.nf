include { valid_chemistries } from './chemistry'
include { resolve_index     } from './indexes'

def samplesheet_required_columns() {
    ['assay', 'experiment_id', 'historical_number', 'replicate',
     'modality', 'chemistry', 'index_type', 'index',
     'species', 'n_donors', 'adt_file']
}

def valid_assays()     { ['GEX', 'CITE', 'DOGMA', 'ATAC', 'Multiome', 'ASAP', 'Flex'] }
def valid_modalities() { ['GEX', 'ATAC', 'ADT', 'HTO', 'VDJ-T', 'VDJ-B', 'CRISPR'] }
def valid_index_types(){ ['SI', 'DI', 'NA'] }

def read_samplesheet_rows(String path) {
    def lines = new File(path).readLines().findAll { line -> !is_blank_row(line) }
    if (lines.isEmpty()) return []
    def headers = lines[0].split(',').collect { h -> h.trim() }
    lines.tail().collect { line ->
        def vals = line.split(',', -1).collect { v -> v.trim() }
        [headers, vals].transpose().collectEntries()
    }
}

def is_blank_row(String line) {
    line.trim().isEmpty() || line.split(',', -1).every { v -> v.trim().isEmpty() }
}

def preflight_samplesheet(String path) {
    def lines = new File(path).readLines()
    if (lines.isEmpty()) error "samplesheet is empty: ${path}"

    def headers = lines[0].split(',').collect { h -> h.trim() }
    def missing = samplesheet_required_columns() - headers
    if (missing) error "samplesheet missing columns: ${missing.join(', ')}"

    def valid_chem = valid_chemistries()

    lines.tail().eachWithIndex { line, i ->
        if (is_blank_row(line)) return
        def vals = line.split(',', -1).collect { v -> v.trim() }
        if (vals.size() != headers.size())
            error "samplesheet row ${i + 2}: expected ${headers.size()} fields, got ${vals.size()}"
        def row = [headers, vals].transpose().collectEntries()

        if (!valid_assays().contains(row.assay))
            error "row ${i + 2}: unknown assay '${row.assay}'. Valid: ${valid_assays().join(', ')}"
        if (!valid_modalities().contains(row.modality))
            error "row ${i + 2}: unknown modality '${row.modality}'. Valid: ${valid_modalities().join(', ')}"
        if (!['human', 'mouse'].any { s -> s.equalsIgnoreCase(row.species) })
            error "row ${i + 2}: unknown species '${row.species}'. Valid: human, mouse"
        if (!valid_index_types().contains(row.index_type))
            error "row ${i + 2}: unknown index_type '${row.index_type}'. Valid: SI, DI, NA"

        if (!valid_chem.contains(row.chemistry))
            error "row ${i + 2}: unrecognised chemistry '${row.chemistry}'. Valid: ${valid_chem.join(', ')}"
    }
}

def resolve_adt_csv(String adt_file, String ss_path, adt_files_dir) {
    if (!adt_file) return null

    def ss_dir     = new File(ss_path).parentFile
    def local_csv  = new File("${ss_dir}/adt_files/${adt_file}.csv")
    def parent_csv = new File("${ss_dir.parentFile}/adt_files/${adt_file}.csv")

    if (local_csv.exists())  return local_csv.canonicalPath
    if (parent_csv.exists()) return parent_csv.canonicalPath
    if (adt_files_dir)       return file("${adt_files_dir}/${adt_file}.csv").toAbsolutePath().toString()

    log.warn "'${adt_file}.csv' not in ${ss_dir}/adt_files/ and --adt_files_dir unset; its library will fail at cellranger multi"
    return null
}

// read_samplesheet_rows has already trimmed every value.
def parse_row(row, Map si_indexes, String ss_path, adt_files_dir) {
    def library_id = "${row.assay}_${row.experiment_id}_exp${row.historical_number}_lib${row.replicate}".toString()
    def n_donors   = row.n_donors.toUpperCase() in ['NA', ''] ? 1 : row.n_donors.toInteger()
    def adt_file   = row.adt_file.toUpperCase() in ['NA', ''] ? null : row.adt_file

    [
        id:                "${library_id}_${row.modality}".toString(),
        library_id:        library_id,
        assay:             row.assay,
        experiment_id:     row.experiment_id,
        historical_number: row.historical_number,
        replicate:         row.replicate,
        modality:          row.modality,
        chemistry:         row.chemistry,
        index_type:        row.index_type,
        index:             row.index,
        index_seqs:        resolve_index(row.index, si_indexes),
        species:           row.species.toLowerCase(),
        n_donors:          n_donors,
        adt_file:          adt_file,
        adt_csv_path:      resolve_adt_csv(adt_file, ss_path, adt_files_dir)
    ]
}

def parse_samplesheet(String ss_path, Map si_indexes, String run_name, adt_files_dir) {
    read_samplesheet_rows(ss_path).collect { row ->
        def meta = parse_row(row, si_indexes, ss_path, adt_files_dir)
        meta.run_name = run_name
        meta
    }
}

def find_meta_conflicts(rows) {
    def fields = ['assay', 'modality', 'chemistry', 'species', 'index_type', 'index']
    rows.groupBy { r -> r.id }.findResults { id, group ->
        def differing = fields.findAll { f ->
            group.collect { r -> r[f] }.unique().size() > 1
        }
        differing
            ? "library '${id}' disagrees across samplesheets on: " +
              differing.collect { f -> "${f}=${group.collect { r -> r[f] }.unique().join(' vs ')}" }.join(', ')
            : null
    }
}
