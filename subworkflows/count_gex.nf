include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas, config_header, adt_csv,
                       //  gex_fastqs, adt_fastqs, hto_fastqs,
                       //  vdj_t_fastqs, vdj_b_fastqs, crispr_fastqs]
    main:
        CELLRANGER_MULTI(ch_libraries)
    emit:
        outs     = CELLRANGER_MULTI.out.outs
}
