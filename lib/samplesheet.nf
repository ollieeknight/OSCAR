// ─── Samplesheet parsing and validation ──────────────────────────────────────

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

// Read a samplesheet into a list of [header: value] maps, skipping blank lines.
// Shared by validation and parsing so both see identical rows.
def read_samplesheet_rows(String path) {
    def lines = new File(path).readLines().findAll { line -> !is_blank_row(line) }
    if (lines.isEmpty()) return []
    def headers = lines[0].split(',').collect { h -> h.trim() }
    lines.tail().collect { line ->
        def vals = line.split(',', -1).collect { v -> v.trim() }
        [headers, vals].transpose().collectEntries()
    }
}

// True for empty lines and for comma-only spacer rows (',,,,,,,,,,'), which are
// commonly used to visually separate experiments in a samplesheet.
def is_blank_row(String line) {
    line.trim().isEmpty() || line.split(',', -1).every { v -> v.trim().isEmpty() }
}

def preflight_samplesheet(String path) {
    def lines = new File(path).readLines()
    if (lines.isEmpty()) error "ERROR: samplesheet is empty: ${path}"

    def headers = lines[0].split(',').collect { h -> h.trim() }
    def missing = samplesheet_required_columns() - headers
    if (missing) error "ERROR: samplesheet missing columns: ${missing.join(', ')}"

    def valid_chem = valid_chemistries()

    lines.tail().eachWithIndex { line, i ->
        if (is_blank_row(line)) return
        def vals = line.split(',', -1).collect { v -> v.trim() }
        if (vals.size() != headers.size())
            error "ERROR: samplesheet row ${i + 2}: expected ${headers.size()} fields, got ${vals.size()}"
        def row = [headers, vals].transpose().collectEntries()

        if (!valid_assays().any { a -> a.equalsIgnoreCase(row.assay) })
            error "ERROR: row ${i + 2}: unknown assay '${row.assay}'. Valid: ${valid_assays().join(', ')}"
        if (!valid_modalities().contains(row.modality))
            error "ERROR: row ${i + 2}: unknown modality '${row.modality}'. Valid: ${valid_modalities().join(', ')}"
        if (!['human', 'mouse'].any { s -> s.equalsIgnoreCase(row.species) })
            error "ERROR: row ${i + 2}: unknown species '${row.species}'. Valid: human, mouse"
        if (!valid_index_types().contains(row.index_type?.trim()))
            error "ERROR: row ${i + 2}: unknown index_type '${row.index_type}'. Valid: SI, DI, NA"

        // Chemistry must match a registered chemistry exactly — the registry in
        // lib/chemistry.nf is the single source of truth.
        def chem = row.chemistry?.trim() ?: ''
        if (!valid_chem.contains(chem))
            error "ERROR: row ${i + 2}: unrecognised chemistry '${row.chemistry}'. Valid: ${valid_chem.join(', ')}"
    }
}

// Resolve the ADT feature-barcode CSV for one row.
// Local-first, centralized fallback:
//   1. {samplesheet_dir}/adt_files/{adt_file}.csv       (co-located with the run)
//   2. {samplesheet_dir}/../adt_files/{adt_file}.csv    (shared across runs)
//   3. {adt_files_dir}/{adt_file}.csv                   (centralized, optional)
def resolve_adt_csv(String adt_file, String ss_path, adt_files_dir) {
    if (!adt_file) return null

    def ss_dir     = new File(ss_path).parentFile
    def local_csv  = new File("${ss_dir}/adt_files/${adt_file}.csv")
    def parent_csv = new File("${ss_dir.parentFile}/adt_files/${adt_file}.csv")

    if (local_csv.exists())  return local_csv.canonicalPath
    if (parent_csv.exists()) return parent_csv.canonicalPath
    if (adt_files_dir)       return file("${adt_files_dir}/${adt_file}.csv").toAbsolutePath().toString()

    log.warn "WARNING: ADT file '${adt_file}.csv' not found at '${ss_dir}/adt_files/' and --adt_files_dir is not set. " +
             "This library will FAIL at cellranger multi (missing [feature] reference). " +
             "Pass --adt_files_dir or place the CSV at '${ss_dir}/adt_files/${adt_file}.csv'."
    return null
}

// Turn one samplesheet row into the meta map used throughout the pipeline.
def parse_row(row, Map si_indexes, String ss_path, adt_files_dir) {
    def n_donors = (row.n_donors == null || row.n_donors.trim() in ['NA', '', 'na']) \
        ? 1 : row.n_donors.trim().toInteger()
    def index    = row.index.trim()
    def raw_adt  = row.adt_file?.trim()
    def adt_file = (raw_adt == null || raw_adt.toUpperCase() == 'NA' || raw_adt == '') ? null : raw_adt

    [
        id:                "${row.assay}_${row.experiment_id}_exp${row.historical_number}_lib${row.replicate}_${row.modality}",
        library_id:        "${row.assay}_${row.experiment_id}_exp${row.historical_number}_lib${row.replicate}",
        assay:             row.assay.trim(),
        experiment_id:     row.experiment_id.trim(),
        historical_number: row.historical_number.trim(),
        replicate:         row.replicate.trim(),
        modality:          row.modality.trim(),
        chemistry:         row.chemistry.trim(),
        index_type:        row.index_type.trim(),
        index:             index,
        index_seqs:        resolve_index(index, si_indexes),
        species:           row.species.trim().toLowerCase(),
        n_donors:          n_donors,
        adt_file:          adt_file,
        adt_csv_path:      resolve_adt_csv(adt_file, ss_path, adt_files_dir)
    ]
}

// Parse a whole samplesheet into meta maps, tagging each with run_name.
def parse_samplesheet(String ss_path, Map si_indexes, String run_name, adt_files_dir) {
    read_samplesheet_rows(ss_path).collect { row ->
        def meta = parse_row(row, si_indexes, ss_path, adt_files_dir)
        meta.run_name = run_name
        meta
    }
}

// Rows sharing a meta.id that disagree on a library-defining field.
//
// With --extra_bcl_dirs the same library is listed in every flowcell's
// samplesheet, and only the first occurrence survives the meta.id dedup in
// main.nf. A typo in a later samplesheet would otherwise be ignored silently
// and the library counted under the wrong chemistry, species or reference.
// Returns a list of human-readable messages, empty when the rows agree.
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
