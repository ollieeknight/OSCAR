// family groups chemistries that share a demux mask; whitelist is the 10x cell barcode list;
// simpleaf is the --chemistry name; flex_bc and flex_wl are the cyto Flex references.
def chemistry_registry() {
    [
        'SC3Pv2':         [family: 'SC3Pv2',  whitelist: '737K-august-2016.txt',         simpleaf: '10xv2'],
        'SC3Pv3':         [family: 'SC3Pv3',  whitelist: '3M-february-2018_TRU.txt.gz',  simpleaf: '10xv3'],
        'SC3Pv4':         [family: 'SC3Pv4',  whitelist: '3M-3pgex-may-2023_TRU.txt.gz', simpleaf: '10xv4-3p'],
        'SC5P':           [family: 'SC5P',    whitelist: '3M-5pgex-jan-2023.txt.gz',     simpleaf: '10xv2-5p'],
        'SC5P-R2':        [family: 'SC5P',    whitelist: '3M-5pgex-jan-2023.txt.gz',     simpleaf: '10xv2-5p'],
        'SC5Pv3':         [family: 'SC5Pv3',  whitelist: '3M-5pgex-jan-2023.txt.gz',     simpleaf: '10xv3-5p'],
        'ARCv1':          [family: 'ARCv1',   whitelist: '737K-arc-v1.txt.gz',           simpleaf: '10xv3'],
        'ARC-v1':         [family: 'ARCv1',   whitelist: '737K-arc-v1.txt.gz',           simpleaf: '10xv3'],
        'Flex-v2-R1':     [family: 'Flex-v2', whitelist: '737K-flex-v2.txt.gz',          simpleaf: '10x-flexv2-gex-3p',
                           flex_bc: 'flex-v2-384.txt', flex_wl: '737K-fixed-rna-profiling.txt.gz'],
        'Flex-v2-RNA-R2': [family: 'Flex-v2', whitelist: '737K-flex-v2.txt.gz',          simpleaf: '10x-flexv2-gex-3p',
                           flex_bc: 'flex-v2-384.txt', flex_wl: '737K-fixed-rna-profiling.txt.gz'],
        'ATAC':           [family: 'ATAC'],
        'NA':             [family: 'NA'],
    ]
}

def valid_chemistries() {
    chemistry_registry().keySet() as List
}

def chemistry_field(String chemistry, String field) {
    def info = chemistry_registry()[chemistry]
    if (!info)
        error "unknown chemistry '${chemistry}'. Valid: ${valid_chemistries().join(', ')}"
    if (!info[field])
        error "chemistry '${chemistry}' has no ${field}"
    info[field]
}

// Flex reads are probe-based, so there are no unspliced reads to count.
def velocity_chemistry(String chemistry) {
    def info = chemistry_registry()[chemistry]
    info?.family == 'Flex-v2' ? null : info?.simpleaf
}
