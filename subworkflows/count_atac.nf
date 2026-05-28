include { CELLRANGER_ATAC } from '../modules/count_atac'

workflow COUNT_ATAC {
    take:
        ch_atac   // [meta, fastq_dir] — one row per ATAC modality

    main:
        CELLRANGER_ATAC(ch_atac)

    emit:
        outs     = CELLRANGER_ATAC.out.outs     // [meta, outs_dir]
}
