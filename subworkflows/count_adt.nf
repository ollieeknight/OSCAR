include { FEATUREMAP      } from '../modules/count_adt'
include { KALLISTO_INDEX  } from '../modules/count_adt'
include { ASAP_TO_KITE    } from '../modules/count_adt'
include { KALLISTO_BUS    } from '../modules/count_adt'
include { BUSTOOLS_CORRECT } from '../modules/count_adt'
include { BUSTOOLS_SORT   } from '../modules/count_adt'
include { BUSTOOLS_COUNT  } from '../modules/count_adt'

// ASAP-only ADT counting pipeline.
// Channel dependency replaces the SLURM --dependency=afterok pattern in the bash pipeline.

workflow COUNT_ADT {
    take:
        // [atac_meta, atac_outs_dir, adt_meta, adt_fastqs]
        // atac_outs must be complete before this subworkflow runs
        ch_asap_adt

    main:
        // Resolve ADT CSV from adt_meta.adt_csv_path (set in parse_row via --adt_files_dir)
        ch_adt_csv = ch_asap_adt
            .map { atac_meta, atac_outs, adt_meta, adt_fastqs ->
                if (!adt_meta.adt_csv_path) {
                    error "Library '${adt_meta.library_id}' has ASAP ADT/HTO modalities but no feature barcode CSV was resolved. " +
                          "Check that 'adt_file' is set in the samplesheet and either place " +
                          "{samplesheet_dir}/adt_files/{adt_file}.csv or pass --adt_files_dir."
                }
                def adt_csv = file(adt_meta.adt_csv_path)
                [ adt_meta, adt_csv, adt_fastqs ]
            }

        FEATUREMAP(ch_adt_csv.map { meta, adt_csv, fqs -> [ meta, adt_csv ] })

        KALLISTO_INDEX(FEATUREMAP.out.index_files)

        // ATAC whitelist is required (set via params.atac_whitelist in nextflow.config)
        ch_whitelist = Channel.value(file(params.atac_whitelist))

        ASAP_TO_KITE(
            ch_adt_csv.map { meta, adt_csv, fqs -> [ meta, fqs ] }
        )

        // Join index + converted FASTQs by meta
        ch_bus_input = KALLISTO_INDEX.out.index
            .join(ASAP_TO_KITE.out.converted_fastqs, by: 0)
            .map { meta, t2g, idx, converted -> [ meta, t2g, idx, converted ] }

        KALLISTO_BUS(ch_bus_input)

        BUSTOOLS_CORRECT(KALLISTO_BUS.out.bus, ch_whitelist)

        BUSTOOLS_SORT(BUSTOOLS_CORRECT.out.corrected)

        BUSTOOLS_COUNT(BUSTOOLS_SORT.out.sorted)

    emit:
        counts   = BUSTOOLS_COUNT.out.counts   // [meta, count_files]
}
