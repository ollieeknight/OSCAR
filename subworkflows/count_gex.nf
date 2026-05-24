include { MULTI_CONFIG }      from '../modules/count_gex'
include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas_list, fastq_dirs, adt_csv]
    main:
        // ← Materialise metas to ArrayList before passing to process
        ch_libraries
            .map { lid, metas, dirs, adt_csv ->
                [lid, metas as ArrayList, dirs as ArrayList, adt_csv]
            }
            .set { ch_libraries_safe }

        MULTI_CONFIG(ch_libraries_safe)
        CELLRANGER_MULTI(MULTI_CONFIG.out.config)
    emit:
        outs     = CELLRANGER_MULTI.out.outs
        versions = CELLRANGER_MULTI.out.versions
}