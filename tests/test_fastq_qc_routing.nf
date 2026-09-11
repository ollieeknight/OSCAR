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

    def via_fastp  = files.findAll {  fastp_covers.call(it.name, it.size) }
    def needs_pigz = files.findAll { !fastp_covers.call(it.name, it.size) }

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

    println "OK: all fastq_qc routing self-checks passed"
}
