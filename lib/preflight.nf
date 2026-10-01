include { preflight_samplesheet; read_samplesheet_rows } from './samplesheet'

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

def resolve_run_from() {
    def rf = (params.run_from ?: 'bcl').toString().toLowerCase()
    if (!['bcl', 'fastq', 'cellranger'].contains(rf))
        error "--run_from must be 'bcl', 'fastq' or 'cellranger' (got: '${params.run_from}')"
    return rf
}

def resolve_extras() {
    def extras = (params.extras ?: '').toString().tokenize(',')*.trim()*.toLowerCase().findAll { e -> e }
    def unknown = extras - ['velocity', 'viral']
    if (unknown)
        error "unknown --extras: ${unknown.join(', ')} (valid: velocity, viral)"
    return extras
}

def preflight(Map paths, String run_from, List extras) {
    check_inputs(paths, run_from, extras)
    ([paths.samplesheet] + paths.extra_samplesheets).each { ss ->
        preflight_samplesheet(ss)
        check_flex(ss)
    }
}

def check_inputs(Map paths, String run_from, List extras) {
    if (!paths.samplesheet) error "--samplesheet is required"
    if (!file(paths.samplesheet).exists()) error "samplesheet not found: ${paths.samplesheet}"

    if (run_from == 'bcl') {
        if (!paths.bcl_dir) error "--bcl_dir is required when --run_from bcl"
        if (!file(paths.bcl_dir).exists()) error "bcl_dir not found: ${paths.bcl_dir}"
    }
    if (run_from == 'fastq') {
        if (!paths.fastq_dir) error "--fastq_dir is required when --run_from fastq"
        paths.fastq_dir.each { d ->
            if (!file(d).exists()) error "fastq_dir not found: ${d}"
        }
    }
    if (run_from == 'cellranger' && !paths.outs_dir)
        error "--outs_dir is required when --run_from cellranger"
    if (run_from == 'bcl' && paths.extra_bcl_dirs && paths.extra_samplesheets
            && paths.extra_samplesheets.size() != paths.extra_bcl_dirs.size())
        error "--extra_samplesheets count (${paths.extra_samplesheets.size()}) must match --extra_bcl_dirs (${paths.extra_bcl_dirs.size()})"

    if ('viral' in extras) {
        if (!params.viral_piscem_index || !file(params.viral_piscem_index).exists())
            error "--extras viral requires viral_piscem_index: ${params.viral_piscem_index}"
        if (!params.viral_t2g || !file(params.viral_t2g).exists())
            error "--extras viral requires viral_t2g: ${params.viral_t2g}"
        if (!(params.bamtofastq_bin ==~ /^https?:.*/) && !file(params.bamtofastq_bin).exists())
            error "--extras viral requires bamtofastq_bin: ${params.bamtofastq_bin}"
    }
    if ('velocity' in extras) {
        def spliceu = [human: params.spliceu_index_human, mouse: params.spliceu_index_mouse]
        def missing = spliceu.findAll { _sp, idx -> !idx || !file(idx).exists() }
        if (missing.size() == spliceu.size())
            error "--extras velocity requires a spliceu index: ${spliceu.values().join(', ')}"
        missing.each { sp, idx -> log.warn "--extras velocity: no ${sp} spliceu index (${idx}), ${sp} libraries will fail" }
    }
}

def check_flex(String ss_path) {
    if (!read_samplesheet_rows(ss_path).any { row -> row.assay == 'Flex' }) return

    if (!(params.flex_backend in ['cellranger', 'cyto', 'both']))
        error "--flex_backend must be 'cellranger', 'cyto', or 'both' (got '${params.flex_backend}')"

    def uses_cr   = params.flex_backend in ['cellranger', 'both']
    def uses_cyto = params.flex_backend in ['cyto', 'both']

    if (uses_cr) {
        if (!params.flex_probe_set)
            error "--flex_backend ${params.flex_backend} needs --flex_probe_set (the standard 10x probe set CSV)"
        if (!file(params.flex_probe_set).exists())
            error "--flex_probe_set not found: ${params.flex_probe_set}"
    }
    if (params.flex_probe_set_custom && !file(params.flex_probe_set_custom).exists())
        error "--flex_probe_set_custom not found: ${params.flex_probe_set_custom}"
    if (params.flex_samples_file && !file(params.flex_samples_file).exists())
        error "--flex_samples_file not found: ${params.flex_samples_file}"

    if (uses_cyto) {
        if (!params.flex_probe_set && !params.flex_probe_set_custom)
            error "--flex_backend cyto needs --flex_probe_set or --flex_probe_set_custom"
        if (!params.container_cyto || !file(params.container_cyto).exists())
            error "cyto container not found at '${params.container_cyto}'. " +
                  "Build it: apptainer build ${params.container_cyto} OSCAR/containers/cyto.def"
    }
}
