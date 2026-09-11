// Self-checks for merged-group samplesheet generation (mixed 8bp SI + 10bp DI
// in ONE bcl-convert call via a per-sample OverrideCycles column).
//
// Exercises the Groovy half of GENERATE_SAMPLESHEET — mask selection per row
// and the global-vs-per-sample decision — plus the grouping/spec-building the
// DEMUX subworkflow does. The shell half (RunInfo.xml wildcard expansion) is
// covered by tests/test_samplesheet_shell.sh.
//
// Run: nextflow run tests/test_samplesheet_gen.nf -w /tmp/nfwork
nextflow.enable.dsl = 2

include { load_si_indexes; resolve_index } from '../lib/indexes'
include { get_override_cycles; apply_si_on_di_correction } from '../modules/demux'

workflow {
    def repo = file("${projectDir}").parent.toString()
    def si   = load_si_indexes(repo, 'novaseq_x')

    // ── apply_si_on_di_correction: index1_cycles is honoured ──────────────
    // Previously hardcoded to 10, which silently produced a short mask on any
    // flow cell whose Index1 read is not 10 cycles.
    assert apply_si_on_di_correction('Y28N*;I10N*;I10N*;Y90N*', 8) == 'Y28N*;I8N2;N*;Y90N*'
    assert apply_si_on_di_correction('Y28N*;I10N*;I10N*;Y90N*', 8, 8)  == 'Y28N*;I8;N*;Y90N*'
    assert apply_si_on_di_correction('Y28N*;I10N*;I10N*;Y90N*', 8, 12) == 'Y28N*;I8N4;N*;Y90N*'
    println "OK: apply_si_on_di_correction honours actual Index1 cycle count"

    // ── the real R462 mix: SC3Pv4 CITE, DI GEX + SI ADT/HTO ───────────────
    // Mirrors the metadata.csv rows: 10bp SI-TT GEX alongside 8bp TruSeq ADT/HTO.
    def rows = [
        [id: 'CITE_A_GEX', modality: 'GEX', index_type: 'DI', index: 'SI-TT-A5'],
        [id: 'CITE_B_GEX', modality: 'GEX', index_type: 'DI', index: 'SI-TT-A6'],
        [id: 'CITE_A_ADT', modality: 'ADT', index_type: 'SI', index: 'ATTCAGAA'],
        [id: 'CITE_A_HTO', modality: 'HTO', index_type: 'SI', index: 'TTCTGAAT'],
    ]
    def metas = rows.collect { r ->
        r + [assay: 'CITE', chemistry: 'SC3Pv4', index_seqs: resolve_index(r.index, si)]
    }

    // Sanity: the mix this whole feature exists for is actually present.
    assert metas.find { it.id == 'CITE_A_GEX' }.index_seqs.is_dual
    assert metas.find { it.id == 'CITE_A_GEX' }.index_seqs.rows[0].i7.length() == 10
    assert !metas.find { it.id == 'CITE_A_ADT' }.index_seqs.is_dual
    assert metas.find { it.id == 'CITE_A_ADT' }.index_seqs.rows[0].i7.length() == 8

    // Same spec-building as subworkflows/demux.nf.
    def is_dual = metas.any { m -> m.index_seqs.is_dual }
    assert is_dual : "group with a DI member must carry an Index2 column"

    def sample_specs = metas.collectMany { m ->
        m.index_seqs.rows.collect { row ->
            [ id: m.id, i7: row.i7, i5: row.get('i5', ''), index_len: row.i7.length(),
              is_dual: m.index_seqs.is_dual, assay: m.assay, chemistry: m.chemistry,
              index_type: m.index_type, modality: m.modality ]
        }
    }

    // Same per-row mask resolution as GENERATE_SAMPLESHEET.
    def specs = sample_specs.collect { sp ->
        def fake = [is_dual: sp.is_dual, rows: [[i7: sp.i7]]]
        def oc4 = get_override_cycles(sp.assay, sp.chemistry, sp.index_type, sp.modality, 4, fake, sp.index_len)
        [id: sp.id, i7: sp.i7, i5: sp.i5, oc4: oc4]
    }
    def by_id = specs.collectEntries { [(it.id): it] }

    // The DI rows keep both 10bp index reads.
    assert by_id['CITE_A_GEX'].oc4 == 'Y28N*;I10N*;I10N*;Y90N*' : by_id['CITE_A_GEX'].oc4

    // The 8bp SI rows use 8 index cycles and fully mask i5. This is the
    // correction that makes one shared call legal — under the old single
    // global 10bp mask bcl-convert aborted with
    // "index of length 8 bases, but a length of 10 was expected".
    assert by_id['CITE_A_ADT'].oc4 == 'Y28N*;I8N2;N*;Y90N*' : by_id['CITE_A_ADT'].oc4
    assert by_id['CITE_A_HTO'].oc4 == 'Y28N*;I8N2;N*;Y90N*' : by_id['CITE_A_HTO'].oc4
    println "OK: mixed group yields distinct per-row masks (I10 for DI, I8N2+masked i5 for SI)"

    // Mixed group => per-sample column required.
    def per_sample = specs.collect { it.oc4 }.unique().size() > 1
    assert per_sample : "mixed 8bp/10bp group must use the per-sample OverrideCycles column"
    println "OK: mixed group selects per-sample OverrideCycles column"

    // ── uniform group keeps the global setting ────────────────────────────
    // Illumina: a setting may be global OR per-sample, never both. A group that
    // needs no per-row variation must not grow the column.
    def uni = ['SI-TT-A5', 'SI-TT-A6', 'SI-TT-A7'].collect { code ->
        [assay: 'CITE', chemistry: 'SC3Pv4', index_type: 'DI', modality: 'GEX',
         index_seqs: resolve_index(code, si)]
    }
    def uni_masks = uni.collect { m ->
        get_override_cycles(m.assay, m.chemistry, m.index_type, m.modality, 4,
                            m.index_seqs, m.index_seqs.rows[0].i7.length())
    }
    assert uni_masks.unique().size() == 1
    println "OK: uniform group keeps a single global OverrideCycles setting"

    // ── index1 column length always matches its own mask ──────────────────
    // The invariant bcl-convert enforces, and the one the reverted attempt
    // violated: for every row, declared I-cycles == length of its Index value.
    specs.each { sp ->
        def i1 = sp.oc4.split(';')[1]
        def used = (i1 =~ /I(\d+)/)[0][1].toInteger()
        assert used == sp.i7.length() :
            "row ${sp.id}: mask ${sp.oc4} declares I${used} but index '${sp.i7}' is ${sp.i7.length()}bp"
    }
    println "OK: every row's declared index cycles match its index length"

    // ── FASTQ matching must not let one id steal another's files ──────────
    // Merging every modality into one group puts all of a flowcell's FASTQ in
    // the same directory, so a bare substring test would mis-assign whenever
    // one library id is a prefix of another.
    def fq_names = [
        'CITE_X_exp1_libA_GEX_S1_L001_R1_001.fastq.gz',
        'CITE_X_exp1_libA_GEX_S1_L001_R2_001.fastq.gz',
        'CITE_X_exp1_libA_GEX_2_S2_L001_R1_001.fastq.gz',
        'CITE_X_exp1_libA_ADT_S3_L001_R1_001.fastq.gz',
        'Undetermined_S0_L001_R1_001.fastq.gz',
    ]
    def match_for = { String id ->
        def re = java.util.regex.Pattern.compile(
            "^" + java.util.regex.Pattern.quote(id) + "_S\\d+_")
        fq_names.findAll { re.matcher(it).find() }
    }
    assert match_for.call('CITE_X_exp1_libA_GEX').size() == 2 : match_for.call('CITE_X_exp1_libA_GEX')
    assert match_for.call('CITE_X_exp1_libA_GEX_2') == ['CITE_X_exp1_libA_GEX_2_S2_L001_R1_001.fastq.gz']
    assert match_for.call('CITE_X_exp1_libA_ADT').size() == 1
    assert match_for.call('CITE_X_exp1_libA_GEX').every { !it.startsWith('CITE_X_exp1_libA_GEX_2') }
    assert match_for.call('CITE_X_exp1_libA').isEmpty() : "partial id must not match anything"
    println "OK: FASTQ matching anchored to bcl-convert naming, no prefix theft"

    // ── merged CITE modalities must share one Y-read structure ────────────
    // Merging is only safe when every modality in the group reads the same
    // cycles as data; only the index part may vary per row. If a CITE mask
    // ever gains a differing Y-structure, merging must stop, so assert it.
    def y_struct = { String oc ->
        def p = oc.split(';') as List
        [p[0], p[-1]]   // read1 and read2 masks
    }
    def cite_masks = ['GEX', 'ADT', 'HTO'].collect { mod ->
        get_override_cycles('CITE', 'SC3Pv4', 'SI', mod, 4,
                            [is_dual: false, rows: [[i7: 'ATTCAGAA']]], 8)
    }
    assert cite_masks.collect { y_struct.call(it) }.unique().size() == 1 :
        "CITE modalities no longer share a Y-read structure: ${cite_masks}"

    // ATAC is the counter-example the grouping key must keep separate: its
    // read1 differs from GEX, so it is deliberately excluded from merging.
    def m_gex  = get_override_cycles('Multiome', 'ARC-v1', 'DI', 'GEX',  4, [is_dual: true, rows: [[i7: 'A'*10]]], 10)
    def m_atac = get_override_cycles('Multiome', 'ARC-v1', 'DI', 'ATAC', 4, [is_dual: true, rows: [[i7: 'A'*8]]], 8)
    assert y_struct.call(m_gex) != y_struct.call(m_atac) :
        "Multiome GEX and ATAC unexpectedly share a Y-read structure"
    println "OK: CITE modalities share a Y-read structure; Multiome ATAC differs and stays unmerged"

    println "OK: all samplesheet-generation self-checks passed"
}
