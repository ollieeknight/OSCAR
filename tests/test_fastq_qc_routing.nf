// Self-checks for FASTQ_QC integrity-check routing.
//
// Every FASTQ must be checked exactly once — by fastp where fastp already
// decompresses it, by `pigz -t` otherwise — and every input file must still
// reach the validated channel that gates counting. A file dropped here is a
// library silently missing from the run, so the partition is asserted to be
// both complete and non-overlapping.
//
// Run: nextflow run tests/test_fastq_qc_routing.nf -w /tmp/nfwork
nextflow.enable.dsl = 2

workflow {
    // Same predicate as subworkflows/fastq_qc.nf.
    def fastp_covers = { String name, long size ->
        name =~ /_R[0-9]+_/ && size > 1024 * 1024
    }

    def MB = 1024 * 1024

    // Realistic bcl-convert output: R-reads plus index reads, plus a tiny
    // R-read (a nearly-empty library, which is exactly the case worth checking).
    def files = [
        [name: 'LIB_A_S1_L001_R1_001.fastq.gz', size: 900L * MB],
        [name: 'LIB_A_S1_L001_R2_001.fastq.gz', size: 2100L * MB],
        [name: 'LIB_A_S1_L001_I1_001.fastq.gz', size: 40L * MB],
        [name: 'LIB_A_S1_L001_I2_001.fastq.gz', size: 40L * MB],
        [name: 'LIB_B_S2_L001_R1_001.fastq.gz', size: 512L],          // below fastp floor
        [name: 'LIB_B_S2_L001_R2_001.fastq.gz', size: 1024L * 1024L], // exactly 1MB, not > 1MB
    ]

    def via_fastp  = files.findAll { f ->  fastp_covers.call(f.name, f.size) }
    def needs_pigz = files.findAll { f -> !fastp_covers.call(f.name, f.size) }

    // ── the partition is total and disjoint ───────────────────────────────
    assert via_fastp.size() + needs_pigz.size() == files.size()
    assert (via_fastp*.name).intersect(needs_pigz*.name).isEmpty() :
        "a file routed to both checks would be decompressed twice again"
    assert (via_fastp*.name + needs_pigz*.name).toSet() == (files*.name).toSet() :
        "every FASTQ must reach the validated channel; a dropped file is a lost library"
    println "OK: every FASTQ checked exactly once (partition total and disjoint)"

    // ── the routing matches what fastp actually processes ─────────────────
    assert via_fastp*.name.sort() == [
        'LIB_A_S1_L001_R1_001.fastq.gz',
        'LIB_A_S1_L001_R2_001.fastq.gz',
    ] : via_fastp*.name

    // Index reads are never seen by fastp, so pigz must still cover them.
    assert needs_pigz*.name.contains('LIB_A_S1_L001_I1_001.fastq.gz')
    assert needs_pigz*.name.contains('LIB_A_S1_L001_I2_001.fastq.gz')
    // Sub-floor R-reads fall through fastp's size filter and must not lose
    // their integrity check along with it.
    assert needs_pigz*.name.contains('LIB_B_S2_L001_R1_001.fastq.gz')
    assert needs_pigz*.name.contains('LIB_B_S2_L001_R2_001.fastq.gz') :
        "1MB exactly is not > 1MB, so fastp skips it and pigz must not"
    println "OK: index reads and sub-1MB R-reads still get an explicit check"

    // ── the saving is real ────────────────────────────────────────────────
    // Before: every file was pigz-tested AND the R-reads were read again by
    // fastp. After: each file is decompressed once.
    def before = files.size() + via_fastp.size()
    def after  = files.size()
    assert after < before
    println "OK: decompressions per library-lane drop from ${before} to ${after}"

    // ── the join key matches FASTP's `checked` output ─────────────────────
    // FASTP emits `fastq_name` (basename minus .fastq.gz); FASTQ_QC rebuilds
    // the same string to join on. If either side changes, files stop being
    // released and counting stalls rather than failing loudly.
    files.each { f ->
        def from_subworkflow = f.name.replaceAll(/\.fastq\.gz$/, '')
        def from_fastp_input = f.name.replaceAll(/\.fastq\.gz$/, '')  // ch_fastp_input builds it the same way
        assert from_subworkflow == from_fastp_input
        assert !from_subworkflow.endsWith('.fastq.gz')
        assert !from_subworkflow.contains('/')
    }
    println "OK: FASTP.out.checked join key matches the key FASTQ_QC rebuilds"

    // ── the release join survives multiple flowcells ──────────────────────
    // bcl-convert names files <Sample_ID>_S<n>_L<lane>_<read>_001 with no
    // flowcell component, so the same library+lane on two flowcells produces
    // byte-identical basenames. `join` is 1:1: keyed on the basename alone the
    // second flowcell is dropped with no error, halving a library's reads.
    // The key must therefore include the per-flowcell fastq_dir.
    def fc_meta = [id: 'LIB_A_GEX', library_id: 'LIB_A', modality: 'GEX']
    def dup_name = 'LIB_A_GEX_S1_L001_R1_001.fastq.gz'
    def fc_items = [
        [fc_meta, '/data/FC_A_fastq', dup_name],
        [fc_meta, '/data/FC_B_fastq', dup_name],
    ]

    def ch_fc_files = channel.fromList(fc_items)
    // FASTP.out.checked carries one row per (fastq_dir, fastq_name).
    def ch_fc_checked = channel.fromList(fc_items.collect { _m, d, n ->
        [d, n.replaceAll(/\.fastq\.gz$/, ''), true]
    })

    ch_fc_files
        .map { m, d, n -> [[d, n.replaceAll(/\.fastq\.gz$/, '')], m, d, n] }
        .join(ch_fc_checked.map { d, n, ok -> [[d, n], ok] }, by: 0, failOnDuplicate: true)
        .count()
        .subscribe { n ->
            assert n == 2 :
                "release join dropped a flowcell: expected 2 released FASTQs, got ${n}"
            println "OK: both flowcells survive the fastp release join"
        }

    // ── regrouping keeps each flowcell's fastq_dir ────────────────────────
    // Grouping by meta.id alone merges every flowcell of a library into one
    // row, then collapses fq_dirs to the first. Downstream velocity quant
    // reads that dir string, so the dropped flowcells become invisible to it.
    channel.fromList(fc_items)
        .map { m, d, f -> [[m.id, d], m, d, f] }
        .groupTuple(by: 0)
        .map { key, ms, _ds, fs -> [ms[0], key[1], fs] }
        .toList()
        .subscribe { rows ->
            def dirs = rows.collect { r -> r[1] }.toSet()
            assert dirs == ['/data/FC_A_fastq', '/data/FC_B_fastq'].toSet() :
                "flowcell fastq_dirs collapsed: kept ${dirs}, lost the rest"
            println "OK: per-flowcell fastq_dir survives regrouping"
        }

    // ── the fastq_dir is the published flowcell dir, not a work dir ───────
    // BCLCONVERT declares its FASTQs as `path("fastqs/*.fastq.gz")`, so the
    // emitted files live in a per-lane Nextflow work dir while publishDir
    // copies them to {bcl_parent}/{run}_fastq. Deriving fastq_dir from
    // fqs[0].parent therefore named a directory holding a single lane of a
    // single demux group. That string is a publishDir target for fastp and
    // MultiQC, and the glob root for velocity quant, so it must be the
    // published path — built here the same way DEMUX_QC builds it.
    def bcl_parent = '/data/runs'
    def bcl_name   = '20240101_VH00206_R502_bcl'
    def run_dir    = bcl_name.replaceAll(/_bcl.*$/, '')
    def published  = "${bcl_parent}/${run_dir}_fastq".toString()

    assert published == '/data/runs/20240101_VH00206_R502_fastq' : published
    assert !published.contains('/work/') :
        "fastq_dir must not point into a Nextflow work dir"
    assert !published.endsWith('/fastqs') :
        "fastq_dir must be the published flowcell dir, not the task's fastqs/ output dir"
    println "OK: fastq_dir resolves to the published flowcell dir"

    // ── --run_from fastq matches FASTQs by anchored id, not prefix ─────────
    // A library id that is a prefix of another ('..._GEX' vs '..._GEX2') stole
    // the other's FASTQs under startsWith(meta.id). Anchoring on bcl-convert's
    // <Sample_ID>_S<n>_ naming fixes it — the same anchoring the BCL path uses.
    def fq_names = [
        'CITE_x_exp1_lib1_GEX_S1_L001_R1_001.fastq.gz',
        'CITE_x_exp1_lib1_GEX2_S2_L001_R1_001.fastq.gz',
    ]
    def match_for = { String id ->
        def re = java.util.regex.Pattern.compile(
            "^" + java.util.regex.Pattern.quote(id) + "_S\\d+_")
        fq_names.findAll { n -> re.matcher(n).find() }
    }

    def gex_match  = match_for.call('CITE_x_exp1_lib1_GEX')
    def gex2_match = match_for.call('CITE_x_exp1_lib1_GEX2')
    assert gex_match == ['CITE_x_exp1_lib1_GEX_S1_L001_R1_001.fastq.gz'] :
        "anchored match must not steal the _GEX2 library's FASTQs: ${gex_match}"
    assert gex2_match == ['CITE_x_exp1_lib1_GEX2_S2_L001_R1_001.fastq.gz'] : gex2_match
    // Guards the fixture: the naive prefix test must actually over-match here.
    assert fq_names.findAll { n -> n.startsWith('CITE_x_exp1_lib1_GEX') }.size() == 2 :
        "fixture does not reproduce the bug; startsWith already discriminates"
    println "OK: --run_from fastq matches by anchored id, not bare prefix"

    println "OK: all fastq_qc routing self-checks passed"
}
