include { CELLSNP_LITE; VIREO } from '../modules/genotype'

workflow GENOTYPE {
    take:
        ch_input
        mode

    main:
        CELLSNP_LITE(ch_input, mode)
        VIREO(CELLSNP_LITE.out.vcf, mode)

    emit:
        VIREO.out.donor_ids
}
