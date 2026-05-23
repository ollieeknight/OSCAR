include { GENERATE_SAMPLESHEET } from '../modules/demux'
include { BCL_TO_FASTQ }        from '../modules/demux'
include { FALCO }               from '../modules/demux'

workflow DEMUX {
    take:
        ch_meta_bcl  // [meta, bcl_dir] pairs — each meta pre-bound to its BCL dir

    main:
        // Group by (assay_indextype_chemistry_modality_indexlength, bcl_dir) → one BCL Convert job per group
        // Index length determines OverrideCycles: 8bp (TruSeq) vs 10bp (10x)
        ch_meta_bcl
            .map { meta, bcl_dir ->
                def index_len = meta.index_seqs?.rows[0]?.i7?.length() ?: 0
                def key = "${meta.assay}_${meta.index_type}_${meta.chemistry}_${meta.modality}_${index_len}_${bcl_dir.name}"
                [key, meta, bcl_dir]
            }
            .groupTuple(by: 0)
            .map { key, metas, bcl_dirs -> [key, metas, bcl_dirs[0]] }
            .set { ch_demux_input }

        GENERATE_SAMPLESHEET(ch_demux_input)
        BCL_TO_FASTQ(GENERATE_SAMPLESHEET.out.samplesheet)

        // Explode output back to individual metas by matching meta.id to FASTQ filename.
        // meta.id == BCL Convert Sample_ID == FASTQ filename prefix.
        // When the same library appears across multiple BCL dirs (multi-run), both
        // sets of FASTQs end up as separate [meta, fastqs] items; groupTuple in main.nf
        // merges them before counting.
        // Derive the published fastq dir path from bcl_dir.name (e.g. R463_bcl → R463_fastq).
        // Pass the directory string rather than staged files — prevents filename collisions
        // when the same library is sequenced on multiple flowcells (identical _S1_ naming).
        BCL_TO_FASTQ.out.fastqs
            .flatMap { metas, bcl_name, fq_files ->
                def fqs     = fq_files instanceof List ? fq_files : [fq_files]
                def run     = bcl_name.replaceAll(/_bcl.*$/, '')
                def fq_dir  = "${params.outdir}/${run}_fastq"
                metas.collectMany { meta ->
                    def matched = fqs.findAll { f -> f.name.contains(meta.id) }
                    matched ? [[meta, fq_dir, matched]] : []
                }
            }
            .set { ch_fastqs }

        // Run Falco per R-read FASTQ (R1/R2/R3 only; I1/I2 index reads skipped)
        // Thread run_name (derived from bcl_dir) so each report lands in the correct
        // {run}_fastq/falco/ directory, not a shared params.run_name dir.
        BCL_TO_FASTQ.out.fastqs
            .flatMap { metas, bcl_name, fq_files ->
                def run = bcl_name.replaceAll(/_bcl.*$/, '')
                (fq_files instanceof List ? fq_files : [fq_files])
                    .findAll { f -> f.name =~ /_R[0-9]+_/ }
                    .collect { f -> [run, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_falco_input }

        FALCO(ch_falco_input)

        FALCO.out.report
            .collect()
            .set { ch_falco_reports }

    emit:
        fastqs        = ch_fastqs        // [meta, fastq_dir_string, [fastq_files]]
        falco_reports = ch_falco_reports // collected falco dirs (passed to MULTIQC in main.nf)
        versions      = BCL_TO_FASTQ.out.versions
}
