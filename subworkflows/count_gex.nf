include { CELLRANGER_MULTI } from '../modules/count_gex'
include { CYTO_FLEX; CYTO_RENAME_SAMPLES } from '../modules/count_gex_cyto'
include { FLEX_PROBE_PREPARE; FLEX_SAMPLE_PREPARE
          FLEX_BARCODE_EXTRACT; FLEX_WHITELIST_EXTRACT } from '../modules/flex_probe_convert'
include { build_multi_config_header; build_flex_samples_section } from '../lib/multi_config'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas, refs, adt_csv,
                       //  gex_fastqs, adt_fastqs, hto_fastqs,
                       //  vdj_t_fastqs, vdj_b_fastqs, crispr_fastqs]

    main:
        // ── Split Flex from regular libraries ─────────────────────────────────
        ch_libraries.branch { entry ->
            flex:    entry[1].any { m -> m.assay == 'Flex' }
            regular: true
        }.set { ch_split }

        def needs_cr   = params.flex_backend in ['cellranger', 'both']
        def needs_cyto = params.flex_backend in ['cyto', 'both']

        // ── Probe set preparation ─────────────────────────────────────────────
        // Runs whenever a custom probe set is supplied, regardless of backend:
        // cellranger needs the merged CSV too, otherwise custom probes are
        // silently absent from its output.
        def has_custom_probes = params.flex_probe_set_custom as boolean

        def ch_probe_csv_cr
        if (has_custom_probes) {
            def std_probe  = params.flex_probe_set
                ? file(params.flex_probe_set)        : file('NO_FILE')
            def cust_probe = file(params.flex_probe_set_custom)

            FLEX_PROBE_PREPARE(channel.value([std_probe, cust_probe]))
            ch_probe_csv_cr = FLEX_PROBE_PREPARE.out.probe_csv_cr
        }
        else {
            // No custom probes — cellranger uses the standard 10x CSV directly.
            ch_probe_csv_cr = params.flex_probe_set
                ? channel.value(file(params.flex_probe_set))
                : channel.value(file('NO_FILE'))
        }

        // ── cellranger path ───────────────────────────────────────────────────
        // Non-Flex always runs cellranger. Flex runs it for backend
        // 'cellranger' or 'both'.
        def ch_for_cr = needs_cr
            ? ch_split.regular.mix(ch_split.flex)
            : ch_split.regular

        // Attach the resolved probe set, then build the config header here —
        // inside the subworkflow, where the merged probe CSV is reachable.
        ch_for_cr
            .combine(ch_probe_csv_cr)
            .map { lid, metas, refs, adt_csv,
                   gex_fqs, adt_fqs, hto_fqs, vdj_t, vdj_b, crispr, probe_csv ->
                def probe_set = probe_csv.name == 'NO_FILE'
                    ? null
                    : probe_csv.toAbsolutePath().toString()
                def header = build_multi_config_header(
                    library_id: lid,
                    metas:      metas,
                    refs:       refs,
                    probe_set:  probe_set,
                    adt_csv:    adt_csv
                )
                def meta            = metas.find { m -> m.modality == 'GEX' } ?: metas[0]
                def samples_section = build_flex_samples_section(meta, params.flex_samples_file)

                [lid, metas, header, adt_csv, samples_section,
                 gex_fqs, adt_fqs, hto_fqs, vdj_t, vdj_b, crispr]
            }
            .set { ch_cr_input }

        CELLRANGER_MULTI(ch_cr_input)

        // ── cyto path ─────────────────────────────────────────────────────────
        // Only when flex_backend == 'cyto' or 'both'. Runs in parallel with
        // cellranger (when 'both'); results used for probe-level QC only.
        if (needs_cyto) {

            // cyto needs the TSV-format probe file, which only
            // FLEX_PROBE_PREPARE produces — run it if it has not run already.
            if (!has_custom_probes) {
                def std_only = params.flex_probe_set
                    ? file(params.flex_probe_set) : file('NO_FILE')
                FLEX_PROBE_PREPARE(channel.value([std_only, file('NO_FILE')]))
            }

            // Detect Flex chemistry from the first Flex library to auto-select
            // barcode ref + preset.
            def ch_flex_chem = ch_split.flex
                .map { _lid, metas, _refs, _adt, _gex, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    def ml = []; metas.each { m -> ml << m }
                    (ml.find { m -> m.modality == 'GEX' }?.chemistry ?: 'Flex-v2-R1')
                }
                .first()

            def ch_cyto_preset = ch_flex_chem.map { chem ->
                chem ==~ /Flex-v2.*/ ? 'gex-v2' : 'gex-v1'
            }

            // Auto-extract probe barcode ref + cell barcode whitelist from the
            // cellranger container.
            FLEX_BARCODE_EXTRACT(ch_flex_chem)
            FLEX_WHITELIST_EXTRACT(ch_flex_chem)

            // Probe barcode expansion — only needed for multiplexed Flex
            def has_samples = params.flex_samples_file as boolean
            def ch_cyto_barcodes

            if (has_samples) {
                FLEX_SAMPLE_PREPARE(
                    channel.value(file(params.flex_samples_file)),
                    FLEX_BARCODE_EXTRACT.out.barcodes
                )
                ch_cyto_barcodes = FLEX_SAMPLE_PREPARE.out.cyto_barcodes
            } else {
                ch_cyto_barcodes = channel.value(file('NO_FILE'))
            }

            // Build cyto input: [lid, metas, probe_tsv, barcodes, whitelist, preset, gex_fastqs]
            ch_split.flex
                .map { lid, metas, _refs, _adt, gex_fqs, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    [lid, metas, gex_fqs]
                }
                .combine(FLEX_PROBE_PREPARE.out.probe_tsv_cyto)
                .combine(ch_cyto_barcodes)
                .combine(FLEX_WHITELIST_EXTRACT.out.whitelist)
                .combine(ch_cyto_preset)
                .map { lid, metas, gex_fqs, probe_tsv, barcodes, whitelist, preset ->
                    [lid, metas, probe_tsv, barcodes, whitelist, preset, gex_fqs]
                }
                .set { ch_cyto_input }

            CYTO_FLEX(ch_cyto_input)

            def ch_samples_file = has_samples
                ? channel.value(file(params.flex_samples_file))
                : channel.value(file('NO_FILE'))

            CYTO_RENAME_SAMPLES(
                CYTO_FLEX.out.counts.combine(ch_samples_file)
            )
        }

    emit:
        // QC runs only on cellranger output — cyto output is probe-level
        // comparison only. When flex_backend == 'cyto', CELLRANGER_MULTI may
        // still run for non-Flex libraries, so its output is always emitted.
        outs = CELLRANGER_MULTI.out.outs
}
