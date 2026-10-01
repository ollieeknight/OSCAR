include { CELLRANGER_MULTI } from '../modules/count_gex'
include { FLEX_PROBE_PREPARE; FLEX_SAMPLE_PREPARE; CYTO_FLEX; CYTO_RENAME_SAMPLES } from '../modules/flex'
include { chemistry_field } from '../lib/chemistry'
include { build_multi_config_header; build_flex_samples_section } from '../lib/multi_config'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas, adt_csv, [gex, adt, hto, vdj_t, vdj_b, crispr FASTQs]]

    main:
        ch_split = ch_libraries.branch { _lid, metas, _adt_csv, _fqs ->
            flex:    metas.any { m -> m.assay == 'Flex' }
            regular: true
        }

        def needs_cr   = params.flex_backend in ['cellranger', 'both']
        def needs_cyto = params.flex_backend in ['cyto', 'both']
        def std_probe  = file(params.flex_probe_set ?: 'NO_FILE')

        if (params.flex_probe_set_custom || needs_cyto)
            FLEX_PROBE_PREPARE(channel.value([std_probe, file(params.flex_probe_set_custom ?: 'NO_FILE')]))

        ch_probe_csv_cr = params.flex_probe_set_custom
            ? FLEX_PROBE_PREPARE.out.probe_csv_cr
            : channel.value(std_probe)

        ch_cr_input = (needs_cr ? ch_split.regular.mix(ch_split.flex) : ch_split.regular)
            .combine(ch_probe_csv_cr)
            .map { lid, metas, adt_csv, fqs, probe_csv ->
                def probe_set = probe_csv.name == 'NO_FILE' ? null : probe_csv.toAbsolutePath().toString()
                def meta      = metas.find { m -> m.modality == 'GEX' } ?: metas[0]
                [lid, meta.run_name, build_multi_config_header(lid, metas, probe_set, adt_csv), adt_csv,
                 build_flex_samples_section(meta, params.flex_samples_file)] + fqs
            }

        CELLRANGER_MULTI(ch_cr_input)

        if (needs_cyto) {
            // toSortedList, not first(): arrival order would change the cache key.
            ch_flex_chem = ch_split.flex
                .map { _lid, metas, _adt_csv, _fqs -> metas.find { m -> m.modality == 'GEX' }?.chemistry ?: 'Flex-v2-R1' }
                .toSortedList()
                .flatMap { chems -> chems.take(1) }

            // [probe barcode translation, cell barcode whitelist, cyto preset]
            ch_flex_refs = ch_flex_chem.map { chem ->
                [file("${params.tenx_barcodes_dir}/translation/${chemistry_field(chem, 'flex_bc')}"),
                 file("${params.tenx_barcodes_dir}/${chemistry_field(chem, 'flex_wl')}"),
                 chem ==~ /Flex-v2.*/ ? 'gex-v2' : 'gex-v1']
            }

            ch_cyto_barcodes = channel.value(file('NO_FILE'))
            if (params.flex_samples_file) {
                FLEX_SAMPLE_PREPARE(channel.value(file(params.flex_samples_file)), ch_flex_refs.map { refs -> refs[0] })
                ch_cyto_barcodes = FLEX_SAMPLE_PREPARE.out.cyto_barcodes
            }

            ch_cyto_input = ch_split.flex
                .map { lid, metas, _adt_csv, fqs -> [lid, metas[0].run_name, fqs[0]] }
                .combine(FLEX_PROBE_PREPARE.out.probe_tsv_cyto)
                .combine(ch_cyto_barcodes)
                .combine(ch_flex_refs.map { _barcodes, whitelist, preset -> [whitelist, preset] })

            CYTO_FLEX(ch_cyto_input)
            CYTO_RENAME_SAMPLES(CYTO_FLEX.out.counts.combine(channel.value(file(params.flex_samples_file ?: 'NO_FILE'))))
        }

    emit:
        // metas stay out of CELLRANGER_MULTI's inputs so editing n_donors or adt_file does not recount.
        ch_libraries.map { lid, metas, _adt_csv, _fqs -> [lid, metas] }.join(CELLRANGER_MULTI.out.outs)
}
