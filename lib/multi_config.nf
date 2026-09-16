def build_multi_config_header(Map opts) {
    def metas      = opts.metas
    def library_id = opts.library_id
    def refs       = opts.refs
    def probe_set  = opts.probe_set
    def adt_csv    = opts.adt_csv

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
            error "Library '${library_id}' is Flex but no probe set was supplied to the multi config builder."
        lines << "probe-set,${probe_set}"
    }

    if (has_vdj) lines += ['', '[vdj]', "reference,${refs.vdj}"]

    if (has_adt) {
        if (!adt_csv || adt_csv.name == 'NO_FILE')
            error "Library '${library_id}' has ADT/HTO modalities but no feature barcode CSV was resolved. " +
                  "Check that 'adt_file' is set in the samplesheet and either place " +
                  "{samplesheet_dir}/adt_files/{adt_file}.csv or pass --adt_files_dir."
        if (!adt_csv.exists())
            error "Library '${library_id}': feature barcode CSV not found: ${adt_csv.toAbsolutePath()}. " +
                  "Check the samplesheet's 'adt_file' column and --adt_files_dir."
        lines += ['', '[feature]', "reference,${adt_csv.toAbsolutePath()}"]
    }

    lines.join('\n')
}

def build_flex_samples_section(Map meta, samples_file) {
    if (meta.assay != 'Flex' || !samples_file) return ''
    def sf = file(samples_file)
    if (!sf.exists()) {
        log.warn "WARNING: --flex_samples_file not found: ${samples_file} — [samples] section will be omitted (singleplex only)"
        return ''
    }
    '\n\n[samples]\n' + sf.text.trim()
}
