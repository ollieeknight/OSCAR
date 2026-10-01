include { FEATUREMAP; KALLISTO_INDEX; ASAP_TO_KITE; KALLISTO_BUS
          BUSTOOLS_CORRECT; BUSTOOLS_SORT; BUSTOOLS_COUNT } from '../modules/count_adt'

workflow COUNT_ADT {
    take:
        ch_asap_adt

    main:
        ch_adt_csv = ch_asap_adt
            .map { _atac_meta, _atac_outs, adt_meta, adt_fastqs ->
                if (!adt_meta.adt_csv_path)
                    error "Library '${adt_meta.library_id}' has ADT/HTO but no adt_file CSV; set adt_file and --adt_files_dir"
                [ adt_meta, file(adt_meta.adt_csv_path), adt_fastqs ]
            }

        FEATUREMAP(ch_adt_csv.map { meta, adt_csv, _fqs -> [ meta, adt_csv ] })

        KALLISTO_INDEX(FEATUREMAP.out.index_files)

        ch_whitelist = channel.value(file(params.atac_whitelist))

        ASAP_TO_KITE(
            ch_adt_csv.map { meta, _adt_csv, fqs -> [ meta, fqs ] }
        )

        KALLISTO_BUS(KALLISTO_INDEX.out.index.join(ASAP_TO_KITE.out.converted_fastqs, by: 0))

        BUSTOOLS_CORRECT(KALLISTO_BUS.out.bus, ch_whitelist)

        BUSTOOLS_SORT(BUSTOOLS_CORRECT.out.corrected)

        BUSTOOLS_COUNT(BUSTOOLS_SORT.out.sorted)

    emit:
        BUSTOOLS_COUNT.out.counts
}
