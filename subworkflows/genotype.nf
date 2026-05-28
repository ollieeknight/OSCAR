include { CELLSNP_LITE } from '../modules/qc'
include { VIREO        } from '../modules/qc'

workflow GENOTYPE {
    take:
        ch_input   // [meta, bam, bai, barcodes]
        mode       // val: 'gex' | 'atac'

    main:
        CELLSNP_LITE(ch_input, mode)
        VIREO(CELLSNP_LITE.out.vcf, mode)

    emit:
        donor_ids = VIREO.out.donor_ids   // [meta, donor_ids.tsv]
}
