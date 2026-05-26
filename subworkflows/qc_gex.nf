include { CELLBENDER; SCRUBLET } from '../modules/qc'
include { CELLSNP_LITE } from '../modules/qc'
include { VIREO       } from '../modules/qc'

workflow QC_GEX {
    take:
        ch_gex_outs   // [library_id, metas, outs_dir]

    main:
        // Flatten to [meta, outs_dir] using first meta (carries species, n_donors etc.)
        ch_gex_outs
            .map { library_id, metas, outs -> [ metas[0] + [library_id: library_id], outs ] }
            .set { ch_input }

        // Ambient RNA removal (GPU)
        CELLBENDER(ch_input)

        // GEX Doublet detection with Scrublet
        SCRUBLET(CELLBENDER.out.h5)

        // Donor demultiplexing — only when n_donors > 1 and species == human
        ch_input
            .filter { meta, outs -> meta.n_donors > 1 && meta.species == 'human' }
            .set { ch_multi_donor }

        // Build [meta, bam, bai, barcodes] for cellsnp-lite
        // BAM: per_sample_outs/{library}/count/sample_alignments.bam
        // Barcodes: from cellbender output_cell_barcodes.csv
        ch_snp_input = ch_multi_donor
            .join(CELLBENDER.out.barcodes, by: 0)
            .map { meta, outs, barcodes ->
                def bam = file("${outs}/per_sample_outs/${meta.library_id}/count/sample_alignments.bam")
                def bai = file("${outs}/per_sample_outs/${meta.library_id}/count/sample_alignments.bam.bai")
                [ meta, bam, bai, barcodes ]
            }

        CELLSNP_LITE(ch_snp_input, 'gex')
        VIREO(CELLSNP_LITE.out.vcf, 'gex')

    emit:
        cellbender = CELLBENDER.out.h5       // [meta, h5]
        doublets   = SCRUBLET.out.doublets   // [meta, doublets.csv]
        vireo      = VIREO.out.donor_ids     // [meta, donor_ids.tsv] (empty if n_donors <= 1)
        logs       = Channel.empty()
        versions   = CELLBENDER.out.versions.mix(SCRUBLET.out.versions)
}
