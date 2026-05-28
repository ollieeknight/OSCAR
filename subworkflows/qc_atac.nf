include { AMULET      } from '../modules/qc'
include { MGATK2      } from '../modules/qc'
include { MACS3       } from '../modules/qc'
include { GENOTYPE    } from './genotype'

workflow QC_ATAC {
    take:
        ch_atac_outs   // [meta, outs_dir]

    main:
        // AMULET doublet detection — always
        AMULET(ch_atac_outs)

        // mgatk2 mitochondrial genotyping — always
        MGATK2(ch_atac_outs)

        // MACS3 custom peak calling — always
        MACS3(ch_atac_outs)

        // Donor demultiplexing — only when n_donors > 1 and species == human
        ch_atac_outs
            .filter { meta, outs -> meta.n_donors > 1 && meta.species == 'human' }
            .set { ch_multi_donor }

        // Build [meta, bam, bai, barcodes] for cellsnp-lite
        // BAM: outs/possorted_bam.bam
        // Barcodes: outs/filtered_peak_bc_matrix/barcodes.tsv
        ch_snp_input = ch_multi_donor
            .map { meta, outs ->
                def bam      = file("${outs}/possorted_bam.bam")
                def bai      = file("${outs}/possorted_bam.bam.bai")
                def barcodes = file("${outs}/filtered_peak_bc_matrix/barcodes.tsv")
                [ meta, bam, bai, barcodes ]
            }

        GENOTYPE(ch_snp_input, 'atac')

    emit:
        amulet   = AMULET.out.summary        // [meta, MultipletSummary.txt]
        mgatk    = MGATK2.out.results        // [meta, mgatk2_out/]
        peaks    = MACS3.out.peaks           // [meta, peaks/]
        vireo    = GENOTYPE.out.donor_ids    // [meta, donor_ids.tsv]
}
