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
        'ARC-v1': [
            family:    'ARCv1',
            whitelist: '737K-arc-v1.txt.gz',
            simpleaf:  '10xv3',
            velocity:  '10xv3',
        ],
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

def valid_chemistries() {
    chemistry_registry().keySet() as List
}

def chemistry_info(String chemistry) {
    def info = chemistry_registry()[chemistry]
    if (!info)
        error "Unknown chemistry '${chemistry}'. Valid: ${valid_chemistries().join(', ')}"
    return info
}

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

def get_chemistry_family(String chemistry) {
    return chemistry_info(chemistry).family
}
