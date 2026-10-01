include { CELLSNP_LITE; VIREO } from '../modules/genotype'

workflow GENOTYPE {
    take:
        ch_input   // [meta, n_donors, bam, bai, barcodes]
        mode       // 'gex' or 'atac'

    main:
        ch_pooled = ch_input
            .filter { _meta, n_donors, _bam, _bai, _bc -> n_donors > 1 }
            .branch { meta, _n, _bam, _bai, _bc ->
                human: meta.species == 'human'
                other: true
            }

        ch_pooled.other.subscribe { meta, n_donors, _bam, _bai, _bc ->
            log.warn "'${meta.library_id}' has n_donors=${n_donors} but species='${meta.species}'; genotyping is human-only, skipping"
        }

        // n_donors stays out of cellsnp's inputs so changing it only reruns vireo.
        CELLSNP_LITE(ch_pooled.human.map { meta, _n, bam, bai, bc -> [meta, bam, bai, bc] }, mode)
        VIREO(CELLSNP_LITE.out.vcf.join(ch_pooled.human.map { meta, n_donors, _bam, _bai, _bc -> [meta, n_donors] }), mode)
}
