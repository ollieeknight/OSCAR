def build_multi_config_header(String library_id, List metas, Map refs, probe_set, adt_csv) {
    def meta    = metas.find { m -> m.modality == 'GEX' } ?: metas[0]
    def has_vdj = metas.any { m -> m.modality in ['VDJ-T', 'VDJ-B'] }
    def has_adt = metas.any { m -> m.modality in ['ADT', 'HTO'] }

    def is_flex_v2 = meta.assay == 'Flex' && meta.chemistry ==~ /Flex-v2.*/

    def lines = ['[gene-expression]']
    if (!is_flex_v2) lines << "reference,${refs.gex}"
    lines << 'create-bam,true'

    if (meta.assay in ['DOGMA', 'Multiome']) {
        lines << 'chemistry,ARC-v1'
    } else if (meta.assay == 'Flex' && meta.chemistry) {
        lines << "chemistry,${meta.chemistry}"
        if (!probe_set)
            error "Library '${library_id}' is Flex but has no probe set"
        lines << "probe-set,${probe_set}"
    }

    if (has_vdj) lines += ['', '[vdj]', "reference,${refs.vdj}"]

    if (has_adt) {
        if (!adt_csv || adt_csv.name == 'NO_FILE')
            error "Library '${library_id}' has ADT/HTO but no adt_file CSV; set adt_file and --adt_files_dir"
        if (!adt_csv.exists())
            error "Library '${library_id}': adt_file CSV not found: ${adt_csv.toAbsolutePath()}"
        lines += ['', '[feature]', "reference,${adt_csv.toAbsolutePath()}"]
    }

    lines.join('\n')
}

// preflight_flex has already checked that samples_file exists.
def build_flex_samples_section(Map meta, samples_file) {
    meta.assay == 'Flex' && samples_file ? '\n\n[samples]\n' + file(samples_file).text.trim() : ''
}
