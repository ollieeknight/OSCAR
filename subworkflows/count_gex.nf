include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas_list, fastq_files (flat), adt_csv]

    main:
        CELLRANGER_MULTI(ch_libraries)

    emit:
        outs     = CELLRANGER_MULTI.out.outs       // [library_id, metas, outs_dir]
        versions = CELLRANGER_MULTI.out.versions
}
