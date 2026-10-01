include { AMULET; MGATK2; MACS3 } from '../modules/qc_atac'
include { GENOTYPE } from './genotype'

workflow QC_ATAC {
    take:
        ch_atac_outs   // [meta, outs]

    main:
        // Only the fields QC uses, so other samplesheet edits do not rerun mgatk2.
        ch_input = ch_atac_outs.map { meta, outs -> [meta.subMap(['library_id', 'run_name', 'species']), meta.n_donors, outs] }
        ch_qc    = ch_input.map { meta, _n, outs -> [meta, outs] }

        AMULET(ch_qc)
        MGATK2(ch_qc)
        MACS3(ch_qc)

        GENOTYPE(
            ch_input.map { meta, n_donors, outs ->
                [meta, n_donors, file("${outs}/possorted_bam.bam"), file("${outs}/possorted_bam.bam.bai"),
                 file("${outs}/filtered_peak_bc_matrix/barcodes.tsv")]
            },
            'atac'
        )
}
