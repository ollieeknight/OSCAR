// ─── Chemistry registry ──────────────────────────────────────────────────────
// Single source of truth for what each OSCAR chemistry string means.
//
// Every consumer (samplesheet validation, viral detection, velocity, Flex
// barcode/whitelist selection) reads from CHEMISTRY below rather than keeping
// its own private lookup table. Adding a chemistry is a one-place edit.
//
// Fields (all optional except `family`):
//   family     — coarse grouping used for OverrideCycles key resolution
//   whitelist  — 10x cell barcode whitelist filename (viral detection)
//   simpleaf   — simpleaf chemistry string (viral detection)
//   velocity   — simpleaf chemistry string for spliced/unspliced quant.
//                Absent = velocity not supported (probe-based chemistries).
//   flex_bc    — Flex probe barcode reference filename (cyto backend)
//   flex_wl    — Flex cell barcode whitelist filename (cyto backend)
//
// The OverrideCycles masks themselves stay in modules/demux.nf: they are keyed
// by assay+index_type+modality as well as chemistry, so they are a demux
// concern rather than a property of the chemistry alone.

def chemistry_registry() {
    [
        'SC3Pv2': [
            family:    'SC3Pv2',
            whitelist: '737K-august-2016.txt',
            simpleaf:  '10xv2',
            velocity:  '10xv2',
        ],
        'SC3Pv3': [
            family:    'SC3Pv3',
            whitelist: '3M-february-2018_TRU.txt.gz',
            simpleaf:  '10xv3',
            velocity:  '10xv3',
        ],
        'SC3Pv4': [
            family:    'SC3Pv4',
            whitelist: '3M-3pgex-may-2023_TRU.txt.gz',
            simpleaf:  '10xv4-3p',
            velocity:  '10xv4-3p',
        ],
        'SC5P': [
            family:    'SC5P',
            whitelist: '3M-5pgex-jan-2023.txt.gz',
            simpleaf:  '10xv2-5p',
            velocity:  '10xv2-5p',
        ],
        // 5' v2 kits are written as SC5P-R2 (read 2) in existing samplesheets.
        'SC5P-R2': [
            family:    'SC5P',
            whitelist: '3M-5pgex-jan-2023.txt.gz',
            simpleaf:  '10xv2-5p',
            velocity:  '10xv2-5p',
        ],
        'SC5Pv3': [
            family:    'SC5Pv3',
            whitelist: '3M-5pgex-jan-2023.txt.gz',
            simpleaf:  '10xv3-5p',
            velocity:  '10xv3-5p',
        ],
        'ARCv1': [
            family:    'ARCv1',
            whitelist: '737K-arc-v1.txt.gz',
            simpleaf:  '10xv3',
            velocity:  '10xv3',
        ],
        // 10x writes this chemistry as 'ARC-v1'; both spellings appear in
        // existing samplesheets, so both are registered.
        'ARC-v1': [
            family:    'ARCv1',
            whitelist: '737K-arc-v1.txt.gz',
            simpleaf:  '10xv3',
            velocity:  '10xv3',
        ],
        // Flex (Fixed RNA Profiling) is probe-based: no intronic signal, so no
        // velocity entry. Both R1 and RNA-R2 variants share the v2 barcode set.
        'Flex-v2-R1': [
            family:    'Flex-v2',
            whitelist: '737K-flex-v2.txt.gz',
            simpleaf:  '10x-flexv2-gex-3p',
            flex_bc:   'flex-v2-384.txt',
            flex_wl:   '737K-fixed-rna-profiling.txt.gz',
        ],
        'Flex-v2-RNA-R2': [
            family:    'Flex-v2',
            whitelist: '737K-flex-v2.txt.gz',
            simpleaf:  '10x-flexv2-gex-3p',
            flex_bc:   'flex-v2-384.txt',
            flex_wl:   '737K-fixed-rna-profiling.txt.gz',
        ],
        'ATAC': [
            family:    'ATAC',
        ],
        'NA': [
            family:    'NA',
        ],
    ]
}

// Valid chemistry names for samplesheet validation — derived, never hand-listed.
def valid_chemistries() {
    chemistry_registry().keySet() as List
}

// Look up one chemistry, erroring with the full valid list on a miss.
def chemistry_info(String chemistry) {
    def info = chemistry_registry()[chemistry]
    if (!info)
        error "Unknown chemistry '${chemistry}'. Valid: ${valid_chemistries().join(', ')}"
    return info
}

// ─── Field accessors ─────────────────────────────────────────────────────────
// Each errors on a chemistry that exists but lacks the field, so a half-added
// chemistry fails loudly at the point of use rather than silently doing nothing.

def get_viral_whitelist(String chemistry, barcodes_dir) {
    def fn = chemistry_info(chemistry).whitelist
    if (!fn) error "VIRAL_DETECT: no whitelist registered for chemistry '${chemistry}'"
    return "${barcodes_dir}/${fn}"
}

def get_simpleaf_chemistry(String chemistry) {
    def s = chemistry_info(chemistry).simpleaf
    if (!s) error "VIRAL_DETECT: no simpleaf chemistry registered for '${chemistry}'"
    return s
}

// Returns null for chemistries where velocity is not meaningful (Flex, ATAC,
// NA) — callers filter these out rather than failing.
def get_velocity_chemistry(String chemistry) {
    return chemistry_registry()[chemistry]?.velocity
}

def get_flex_barcode_file(String chemistry) {
    def fn = chemistry_info(chemistry).flex_bc
    if (!fn) error "Flex backend: chemistry '${chemistry}' has no probe barcode reference"
    return fn
}

def get_flex_whitelist_file(String chemistry) {
    def fn = chemistry_info(chemistry).flex_wl
    if (!fn) error "Flex backend: chemistry '${chemistry}' has no cell barcode whitelist"
    return fn
}

// Coarse family, used by demux to resolve OverrideCycles keys.
def get_chemistry_family(String chemistry) {
    return chemistry_info(chemistry).family
}
