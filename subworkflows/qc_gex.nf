include { CELLBENDER; SCRUBLET } from '../modules/qc_gex'
include { GENOTYPE } from './genotype'

workflow QC_GEX {
    take:
        ch_gex_outs

    main:
        ch_gex_outs
            .map { library_id, metas, outs -> [ metas[0] + [library_id: library_id], outs ] }
            .set { ch_input }

        CELLBENDER(ch_input)
        SCRUBLET(CELLBENDER.out.h5)

        ch_input
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species == 'human' }
            .set { ch_multi_donor }

        ch_input
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species != 'human' }
            .subscribe { meta, _outs ->
                log.warn "WARN: '${meta.library_id}' declares n_donors=${meta.n_donors} but species='${meta.species}' — genotyping is human-only, skipping donor demultiplexing"
            }

        ch_snp_input = ch_multi_donor
            .join(CELLBENDER.out.barcodes, by: 0)
            .map { meta, outs, barcodes ->
                def bam = file("${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam")
                def bai = file("${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam.bai")
                [ meta, bam, bai, barcodes ]
            }

        GENOTYPE(ch_snp_input, 'gex')

    emit:
        cellbender = CELLBENDER.out.h5
        barcodes   = CELLBENDER.out.barcodes
        doublets   = SCRUBLET.out.doublets
        vireo      = GENOTYPE.out
}
