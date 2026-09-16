include { CELLRANGER_ATAC } from '../modules/count_atac'

workflow COUNT_ATAC {
    take:
        ch_atac_libraries

    main:
        CELLRANGER_ATAC(ch_atac_libraries)

    emit:
        CELLRANGER_ATAC.out.outs
}
