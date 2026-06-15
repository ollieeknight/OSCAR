include { CELLRANGER_MULTI    } from '../modules/count_gex'
include { CYTO_FLEX           } from '../modules/count_gex_cyto'
include { CYTO_RENAME_SAMPLES } from '../modules/count_gex_cyto'
include { FLEX_PROBE_PREPARE    } from '../modules/flex_probe_convert'
include { FLEX_SAMPLE_PREPARE   } from '../modules/flex_probe_convert'
include { FLEX_BARCODE_EXTRACT  } from '../modules/flex_probe_convert'
include { FLEX_WHITELIST_EXTRACT } from '../modules/flex_probe_convert'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas, config_header, adt_csv, flex_samples_content,
                       //  gex_fastqs, adt_fastqs, hto_fastqs,
                       //  vdj_t_fastqs, vdj_b_fastqs, crispr_fastqs]

    main:
        // ── Split Flex from regular libraries ─────────────────────────────────
        ch_libraries.branch {
            flex:    it[1].any { m -> m.assay == 'Flex' }
            regular: true
        }.set { ch_split }

        // ── cellranger path ───────────────────────────────────────────────────
        // Non-Flex always runs cellranger.
        // Flex runs cellranger when backend is 'cellranger' or 'both'.
        def ch_for_cr = (params.flex_backend in ['cellranger', 'both'])
            ? ch_split.regular.mix(ch_split.flex)
            : ch_split.regular

        CELLRANGER_MULTI(ch_for_cr)

        // ── cyto path ─────────────────────────────────────────────────────────
        // Only when flex_backend == 'cyto' or 'both'. Runs in parallel with
        // cellranger (when 'both'); results used for probe-level QC only.
        if (params.flex_backend in ['cyto', 'both']) {

            // Probe format conversion — runs once regardless of library count
            def std_probe  = params.flex_probe_set
                ? file(params.flex_probe_set)        : file('NO_FILE')
            def cust_probe = params.flex_probe_set_custom
                ? file(params.flex_probe_set_custom) : file('NO_FILE')

            FLEX_PROBE_PREPARE(Channel.value([std_probe, cust_probe]))

            // Detect Flex chemistry version from first Flex library to auto-select barcode ref + preset
            def ch_flex_chem = ch_split.flex
                .map { lid, metas, _cfg, _adt, _flex, _gex, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    def ml = []; metas.each { ml << it }
                    (ml.find { it.modality == 'GEX' }?.chemistry ?: 'Flex-v2-R1')
                }
                .first()

            def ch_cyto_preset = ch_flex_chem.map { chem ->
                chem ==~ /Flex-v2.*/ ? 'gex-v2' : 'gex-v1'
            }

            // Auto-extract probe barcode ref + cell barcode whitelist from cellranger container
            FLEX_BARCODE_EXTRACT(ch_flex_chem)
            FLEX_WHITELIST_EXTRACT(ch_flex_chem)

            // Probe barcode expansion — only needed for multiplexed Flex
            def has_samples = params.flex_samples_file as boolean
            def ch_cyto_barcodes

            if (has_samples) {
                FLEX_SAMPLE_PREPARE(
                    Channel.value(file(params.flex_samples_file)),
                    FLEX_BARCODE_EXTRACT.out.barcodes.first()
                )
                ch_cyto_barcodes = FLEX_SAMPLE_PREPARE.out.cyto_barcodes.first()
            } else {
                ch_cyto_barcodes = Channel.value(file('NO_FILE'))
            }

            // Build cyto input: [lid, metas, probe_tsv, barcodes, whitelist, cyto_preset, gex_fastqs]
            ch_split.flex
                .map { lid, metas, _cfg, _adt, _flex, gex_fqs, _adt_fqs, _hto_fqs, _vdj_t, _vdj_b, _crispr ->
                    [lid, metas, gex_fqs]
                }
                .combine(FLEX_PROBE_PREPARE.out.probe_tsv_cyto)
                .combine(ch_cyto_barcodes)
                .combine(FLEX_WHITELIST_EXTRACT.out.whitelist.first())
                .combine(ch_cyto_preset)
                .map { lid, metas, gex_fqs, probe_tsv, barcodes, whitelist, preset ->
                    [lid, metas, probe_tsv, barcodes, whitelist, preset, gex_fqs]
                }
                .set { ch_cyto_input }

            CYTO_FLEX(ch_cyto_input)

            def ch_samples_file = has_samples
                ? Channel.value(file(params.flex_samples_file))
                : Channel.value(file('NO_FILE'))

            CYTO_RENAME_SAMPLES(
                CYTO_FLEX.out.counts.combine(ch_samples_file)
            )
        }

    emit:
        // QC runs only on cellranger output — cyto output is probe-level comparison only.
        // When flex_backend == 'cyto', CELLRANGER_MULTI never runs; emit empty channel
        // so downstream QC_GEX is simply a no-op rather than a channel resolution error.
        outs = (params.flex_backend in ['cellranger', 'both'])
            ? CELLRANGER_MULTI.out.outs
            : Channel.empty()
}
