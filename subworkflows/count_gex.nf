include { CELLRANGER_MULTI } from '../modules/count_gex'
include { FLEX_PROBE_PREPARE; FLEX_SAMPLE_PREPARE
          FLEX_REFS; CYTO_FLEX; CYTO_RENAME_SAMPLES } from '../modules/flex'
include { build_multi_config_header; build_flex_samples_section } from '../lib/multi_config'

workflow COUNT_GEX {
    take:
        ch_libraries

    main:
        ch_libraries.branch { entry ->
            flex:    entry[1].any { m -> m.assay == 'Flex' }
            regular: true
        }.set { ch_split }

        def needs_cr   = params.flex_backend in ['cellranger', 'both']
        def needs_cyto = params.flex_backend in ['cyto', 'both']
        def std_probe  = file(params.flex_probe_set ?: 'NO_FILE')

        if (params.flex_probe_set_custom || needs_cyto)
            FLEX_PROBE_PREPARE(channel.value([std_probe, file(params.flex_probe_set_custom ?: 'NO_FILE')]))

        def ch_probe_csv_cr = params.flex_probe_set_custom
            ? FLEX_PROBE_PREPARE.out.probe_csv_cr
            : channel.value(std_probe)

        (needs_cr ? ch_split.regular.mix(ch_split.flex) : ch_split.regular)
            .combine(ch_probe_csv_cr)
            .map { lid, metas, refs, adt_csv,
                   gex_fqs, adt_fqs, hto_fqs, vdj_t, vdj_b, crispr, probe_csv ->
                def probe_set = probe_csv.name == 'NO_FILE' ? null : probe_csv.toAbsolutePath().toString()
                def header    = build_multi_config_header(lid, metas, refs, probe_set, adt_csv)
                def meta      = metas.find { m -> m.modality == 'GEX' } ?: metas[0]

                [lid, meta.run_name, header, adt_csv, build_flex_samples_section(meta, params.flex_samples_file),
                 gex_fqs, adt_fqs, hto_fqs, vdj_t, vdj_b, crispr]
            }
            .set { ch_cr_input }

        CELLRANGER_MULTI(ch_cr_input)

        if (needs_cyto) {
            // toSortedList, not first(): arrival order would change the cache key.
            def ch_flex_chem = ch_split.flex
                .map { _lid, metas, _refs, _adt, _gex, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    metas.find { m -> m.modality == 'GEX' }?.chemistry ?: 'Flex-v2-R1'
                }
                .toSortedList()
                .flatMap { chems -> chems.take(1) }

            def ch_cyto_preset = ch_flex_chem.map { chem -> chem ==~ /Flex-v2.*/ ? 'gex-v2' : 'gex-v1' }

            FLEX_REFS(ch_flex_chem)

            def ch_cyto_barcodes = channel.value(file('NO_FILE'))
            if (params.flex_samples_file) {
                FLEX_SAMPLE_PREPARE(channel.value(file(params.flex_samples_file)), FLEX_REFS.out.barcodes)
                ch_cyto_barcodes = FLEX_SAMPLE_PREPARE.out.cyto_barcodes
            }

            ch_split.flex
                .map { lid, metas, _refs, _adt, gex_fqs, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    [lid, metas[0].run_name, gex_fqs]
                }
                .combine(FLEX_PROBE_PREPARE.out.probe_tsv_cyto)
                .combine(ch_cyto_barcodes)
                .combine(FLEX_REFS.out.whitelist)
                .combine(ch_cyto_preset)
                .map { lid, run_name, gex_fqs, probe_tsv, barcodes, whitelist, preset ->
                    [lid, run_name, probe_tsv, barcodes, whitelist, preset, gex_fqs]
                }
                .set { ch_cyto_input }

            CYTO_FLEX(ch_cyto_input)

            CYTO_RENAME_SAMPLES(
                CYTO_FLEX.out.counts.combine(channel.value(file(params.flex_samples_file ?: 'NO_FILE')))
            )
        }

    emit:
        // metas stay out of CELLRANGER_MULTI's inputs so editing n_donors or adt_file does not recount.
        ch_libraries.map { entry -> [entry[0], entry[1]] }.join(CELLRANGER_MULTI.out.outs)
}
