// Self-checks for lib/ helpers. Run:
//   nextflow run tests/test_lib.nf -w /tmp/nfwork
nextflow.enable.dsl = 2

include { load_si_indexes; resolve_index } from '../lib/indexes'
include { chemistry_registry; valid_chemistries; get_velocity_chemistry
          get_simpleaf_chemistry; get_viral_whitelist
          get_flex_barcode_file; get_flex_whitelist_file
          get_chemistry_family } from '../lib/chemistry'
include { preflight_samplesheet; read_samplesheet_rows } from '../lib/samplesheet'
include { to_abs_path; to_abs_list } from '../main'
include { build_multi_config_header; build_flex_samples_section } from '../lib/multi_config'

workflow {
    // ── indexes ───────────────────────────────────────────────────────────
    // projectDir is tests/; load_si_indexes appends /assets/indexes to the
    // path it is given, so hand it the repo root.
    def repo = file("${projectDir}").parent.toString()
    def si = load_si_indexes(repo, 'novaseq_x')
    assert !si.dual.containsKey('index_name') : "TS header row leaked into dual map"
    assert si.dual['SI-TT-A1'].i7 == 'GTAACATGCG'
    assert si.dual['SI-TT-A1'].i5 == 'AGTGTTACCT'   // workflow A (fwd)
    assert si.dual['SI-TS-A1'].i7 == 'AATTTCGGGT'
    assert si.single['SI-GA-A1'].size() == 4
    def si6k = load_si_indexes(repo, 'novaseq6000')
    assert si6k.dual['SI-TT-A1'].i5 == 'AGGTAACACT' // workflow B (RC)

    // resolve_index: kit code, raw sequence, NA
    assert resolve_index('SI-TT-A1', si).is_dual
    assert resolve_index('SI-GA-A1', si).rows.size() == 4
    assert !resolve_index('SI-GA-A1', si).is_dual
    assert resolve_index('ATTACTCG', si).rows[0].i7 == 'ATTACTCG'
    assert resolve_index('NA', si).rows.isEmpty()

    // ── chemistry registry ────────────────────────────────────────────────
    // Every chemistry accepted by the samplesheet must resolve to a family.
    valid_chemistries().each { c ->
        assert get_chemistry_family(c) != null : "no family for chemistry '${c}'"
    }
    // Flex-v2-RNA-R2 was previously absent from the velocity/viral maps.
    assert 'Flex-v2-RNA-R2' in valid_chemistries()
    assert get_simpleaf_chemistry('Flex-v2-RNA-R2') == '10x-flexv2-gex-3p'
    assert get_velocity_chemistry('Flex-v2-RNA-R2') == null   // probe-based: no velocity
    assert get_velocity_chemistry('SC3Pv3') == '10xv3'
    assert get_viral_whitelist('SC3Pv3', '/bc').endsWith('3M-february-2018_TRU.txt.gz')
    assert get_flex_barcode_file('Flex-v2-R1') == 'flex-v2-384.txt'
    assert get_flex_whitelist_file('Flex-v2-RNA-R2') == '737K-fixed-rna-profiling.txt.gz'

    // ── multi config ──────────────────────────────────────────────────────
    def refs = [gex: '/ref/gex', vdj: '/ref/vdj']

    // plain GEX: reference + create-bam, no probe-set, no feature block
    def gex_meta = [modality: 'GEX', assay: 'GEX', chemistry: 'SC3Pv3']
    def h1 = build_multi_config_header(
        library_id: 'L1', metas: [gex_meta], refs: refs, probe_set: null, adt_csv: null)
    assert h1.contains('reference,/ref/gex')
    assert h1.contains('create-bam,true')
    assert !h1.contains('probe-set')
    assert !h1.contains('[feature]')

    // CITE: ADT present → [feature] block with the resolved CSV
    def adt_meta = [modality: 'ADT', assay: 'CITE', chemistry: 'SC3Pv3']
    def h2 = build_multi_config_header(
        library_id: 'L2', metas: [gex_meta, adt_meta], refs: refs,
        probe_set: null, adt_csv: file('/tmp/adt.csv'))
    assert h2.contains('[feature]')
    assert h2.contains('/tmp/adt.csv')

    // Flex v2: no transcriptome reference, probe-set present.
    // This is the custom-probe path: a merged CSV must land in probe-set.
    def flex_meta = [modality: 'GEX', assay: 'Flex', chemistry: 'Flex-v2-R1']
    def h3 = build_multi_config_header(
        library_id: 'L3', metas: [flex_meta], refs: refs,
        probe_set: '/work/merged_probes_cellranger.csv', adt_csv: null)
    assert !h3.contains('reference,/ref/gex') : "Flex v2 must not set a transcriptome reference"
    assert h3.contains('probe-set,/work/merged_probes_cellranger.csv')
    assert h3.contains('chemistry,Flex-v2-R1')

    // VDJ present → [vdj] block
    def vdj_meta = [modality: 'VDJ-T', assay: 'GEX', chemistry: 'SC5P']
    def h4 = build_multi_config_header(
        library_id: 'L4', metas: [gex_meta, vdj_meta], refs: refs,
        probe_set: null, adt_csv: null)
    assert h4.contains('[vdj]')
    assert h4.contains('reference,/ref/vdj')

    // ADT without a resolved CSV must fail loudly, not silently omit [feature]
    def threw = false
    try {
        build_multi_config_header(
            library_id: 'L5', metas: [gex_meta, adt_meta], refs: refs,
            probe_set: null, adt_csv: file('NO_FILE'))
    } catch (Exception e) { threw = true }
    assert threw : "ADT with no feature CSV must error"

    // singleplex Flex → no [samples] section
    assert build_flex_samples_section(flex_meta, null) == ''

    // The shipped example samplesheet must pass validation, and every chemistry
    // it uses must be registered. Comma-only spacer rows must be skipped.
    def ex = "${repo}/assets/example_metadata.csv"
    if (file(ex).exists()) {
        preflight_samplesheet(ex)
        def rows = read_samplesheet_rows(ex)
        assert rows.every { r -> r.chemistry in valid_chemistries() }
        assert rows.every { r -> r.assay?.trim() } : "spacer row leaked into parsed rows"
        assert rows.size() == 19 : "expected 19 real rows, got ${rows.size()}"
    }

    // Path resolution: relative in, absolute out; comma lists become Lists.
    assert to_abs_path(null) == null
    assert to_abs_path('assets').startsWith('/')
    assert to_abs_list(null) == []
    assert to_abs_list('a,b').size() == 2
    assert to_abs_list('a, b').every { pth -> pth.startsWith('/') }

    // The browser tool's dropdowns are generated from these same definitions
    // (assets/generate_tool_schema.py). Assert the generated schema is present
    // and in step, so a chemistry added to the registry without regenerating
    // the schema fails here rather than on someone's samplesheet.
    def schema_js = file("${repo}/docs-src/tools/helpers/schema.js")
    if (schema_js.exists()) {
        def txt = schema_js.text
        valid_chemistries().each { c ->
            assert txt.contains("\"${c}\"") : "schema.js is stale: missing chemistry '${c}'. Run: python3 assets/generate_tool_schema.py"
        }
    }

    println "OK: all lib/ self-checks passed (${valid_chemistries().size()} chemistries, ${si.dual.size()} dual indexes)"
}
