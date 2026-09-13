include { CELLRANGER_ATAC } from '../modules/count_atac'

workflow COUNT_ATAC {
    take:
        ch_atac_libraries   // [meta, [atac_fastqs]]

    main:
        CELLRANGER_ATAC(ch_atac_libraries)

    emit:
        CELLRANGER_ATAC.out.outs   // [meta, outs/]
}
