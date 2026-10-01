include { CELLBENDER; SCRUBLET } from '../modules/qc_gex'
include { GENOTYPE } from './genotype'

workflow QC_GEX {
    take:
        ch_gex_outs   // [library_id, metas, outs]

    main:
        // Only the fields QC uses, so other samplesheet edits do not rerun cellbender.
        ch_input = ch_gex_outs.map { _lid, metas, outs ->
            [metas[0].subMap(['library_id', 'run_name', 'species']), metas[0].n_donors, outs]
        }

        CELLBENDER(ch_input.map { meta, _n, outs -> [meta, outs] })
        SCRUBLET(CELLBENDER.out.h5)

        GENOTYPE(
            ch_input
                .join(CELLBENDER.out.barcodes)
                .map { meta, n_donors, outs, barcodes ->
                    def bam = "${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam"
                    [meta, n_donors, file(bam), file("${bam}.bai"), barcodes]
                },
            'gex'
        )

    emit:
        CELLBENDER.out.barcodes
}
