include { FEATUREMAP; ASAP_TO_KITE; KALLISTO_BUS; BUSTOOLS } from '../modules/count_adt'

workflow COUNT_ADT {
    take:
        ch_adt   // [meta, fastqs]

    main:
        FEATUREMAP(ch_adt.map { meta, _fqs ->
            if (!meta.adt_csv_path)
                error "library '${meta.library_id}' has ADT/HTO but no adt_file CSV; set adt_file and --adt_files_dir"
            [meta, file(meta.adt_csv_path)]
        })
        ASAP_TO_KITE(ch_adt)
        KALLISTO_BUS(FEATUREMAP.out.features.join(ASAP_TO_KITE.out.fastqs))
        BUSTOOLS(KALLISTO_BUS.out.bus, channel.value(file(params.atac_whitelist)))
}
