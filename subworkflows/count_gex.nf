include { MULTI_CONFIG }      from '../modules/count_gex'
include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas_list, fastq_dirs, adt_csv]

    main:
        MULTI_CONFIG(ch_libraries)
        CELLRANGER_MULTI(MULTI_CONFIG.out.config)

    emit:
        outs     = CELLRANGER_MULTI.out.outs       // [library_id, metas, outs_dir]
        versions = CELLRANGER_MULTI.out.versions
}
