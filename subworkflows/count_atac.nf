include { CELLRANGER_ATAC } from '../modules/count_atac'

workflow COUNT_ATAC {
    take:
        ch_atac_libraries   // [meta, fastqs]

    main:
        // Only the fields cellranger-atac uses, so editing n_donors or adt_file does not recount.
        CELLRANGER_ATAC(ch_atac_libraries.map { meta, fqs ->
            [meta.subMap(['id', 'library_id', 'assay', 'species', 'run_name']), fqs]
        })

    emit:
        ch_atac_libraries.map { meta, _fqs -> [meta.library_id, meta] }
            .join(CELLRANGER_ATAC.out.outs.map { meta, outs -> [meta.library_id, outs] })
            .map { _lid, meta, outs -> [meta, outs] }
}
