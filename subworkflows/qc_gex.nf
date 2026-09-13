include { CELLBENDER; SCRUBLET } from '../modules/qc_gex'
include { GENOTYPE } from './genotype'

workflow QC_GEX {
    take:
        ch_gex_outs   // [library_id, metas, outs_dir]

    main:
        // Flatten to [meta, outs_dir] using first meta (carries species, n_donors etc.)
        ch_gex_outs
            .map { library_id, metas, outs -> [ metas[0] + [library_id: library_id], outs ] }
            .set { ch_input }

        CELLBENDER(ch_input)
        SCRUBLET(CELLBENDER.out.h5)

        // Donor demultiplexing — human only: vireo needs a human SNP reference
        // panel, so a mouse library never takes this path regardless of n_donors.
        // Warn when a samplesheet asks for donors that cannot be resolved, rather
        // than dropping the library silently.
        ch_input
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species == 'human' }
            .set { ch_multi_donor }

        ch_input
            .filter { meta, _outs -> meta.n_donors > 1 && meta.species != 'human' }
            .subscribe { meta, _outs ->
                log.warn "WARN: '${meta.library_id}' declares n_donors=${meta.n_donors} but species='${meta.species}' — genotyping is human-only, skipping donor demultiplexing"
            }

        // cellsnp-lite input: the per-sample BAM plus cellbender's cell list.
        ch_snp_input = ch_multi_donor
            .join(CELLBENDER.out.barcodes, by: 0)
            .map { meta, outs, barcodes ->
                def bam = file("${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam")
                def bai = file("${outs}/per_sample_outs/${meta.library_id}/sample_alignments.bam.bai")
                [ meta, bam, bai, barcodes ]
            }

        GENOTYPE(ch_snp_input, 'gex')

    emit:
        cellbender = CELLBENDER.out.h5        // [meta, h5]
        barcodes   = CELLBENDER.out.barcodes  // [meta, output_cell_barcodes.csv]
        doublets   = SCRUBLET.out.doublets    // [meta, doublets.csv]
        vireo      = GENOTYPE.out   // [meta, donor_ids.tsv] (empty if n_donors <= 1)
}
