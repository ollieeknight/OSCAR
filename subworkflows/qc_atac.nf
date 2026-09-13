include { AMULET; MGATK2; MACS3 } from '../modules/qc_atac'
include { GENOTYPE } from './genotype'

workflow QC_ATAC {
    take:
        ch_atac_outs   // [meta, outs_dir]

    main:
        AMULET(ch_atac_outs)
        MGATK2(ch_atac_outs)
        MACS3(ch_atac_outs)

        // Donor demultiplexing — human only: vireo needs a human SNP reference
        // panel, so a mouse library never takes this path regardless of n_donors.
        // Warn when a samplesheet asks for donors that cannot be resolved, rather
        // than dropping the library silently.
        ch_atac_outs
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species == 'human' }
            .set { ch_multi_donor }

        ch_atac_outs
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species != 'human' }
            .subscribe { meta, _outs ->
                log.warn "WARN: '${meta.library_id}' declares n_donors=${meta.n_donors} but species='${meta.species}' — genotyping is human-only, skipping donor demultiplexing"
            }

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
        vireo    = GENOTYPE.out    // [meta, donor_ids.tsv]
}
