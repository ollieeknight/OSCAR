include { CELLBENDER; SCRUBLET } from '../modules/qc_gex'
include { GENOTYPE } from './genotype'

workflow QC_GEX {
    take:
        ch_gex_outs

    main:
        // Only the fields QC uses, so other samplesheet edits do not rerun cellbender.
        ch_gex_outs
            .map { library_id, metas, outs ->
                [[library_id: library_id, run_name: metas[0].run_name, species: metas[0].species], metas[0].n_donors, outs]
            }
            .set { ch_input }

        CELLBENDER(ch_input.map { meta, _n, outs -> [meta, outs] })
        SCRUBLET(CELLBENDER.out.h5)

        ch_input
            .join(CELLBENDER.out.barcodes)
            .map { meta, n_donors, outs, barcodes ->
                def bam = "${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam"
                [meta, n_donors, file(bam), file("${bam}.bai"), barcodes]
            }
            .set { ch_snp_input }

        GENOTYPE(ch_snp_input, 'gex')

    emit:
        CELLBENDER.out.barcodes
}
