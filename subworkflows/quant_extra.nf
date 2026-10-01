include { VIRAL_DETECT; SIMPLEAF_VELOCITY } from '../modules/quant_extra'
include { chemistry_field; velocity_chemistry } from '../lib/chemistry'

workflow QUANT_EXTRA {
    take:
        ch_gex_outs        // [library_id, metas, outs]
        ch_gex_fastqs      // [meta, fastq_dir, fastqs]
        ch_cellbender_bc   // [meta, barcodes]
        extras

    main:
        if ('viral' in extras) {
            ch_viral_input = ch_gex_outs
                .filter { _lid, metas, _outs -> metas[0].species == 'human' }
                .map { lid, metas, outs ->
                    def chem = metas[0].chemistry
                    [[library_id: lid, run_name: metas[0].run_name],
                     file("${outs}/unassigned_alignments.bam"), file("${outs}/unassigned_alignments.bam.bai"),
                     file("${params.tenx_barcodes_dir}/${chemistry_field(chem, 'whitelist')}"), chemistry_field(chem, 'simpleaf')]
                }

            VIRAL_DETECT(ch_viral_input, file(params.viral_piscem_index), file(params.viral_t2g), file(params.bamtofastq_bin))
        }

        if ('velocity' in extras) {
            ch_velocity_input = ch_gex_fastqs
                .filter { meta, _fastq_dir, _fqs -> meta.modality == 'GEX' && velocity_chemistry(meta.chemistry) }
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, [fastq_dir, fqs.findAll { f -> f.name =~ /_R[12]_/ }]] }
                .groupTuple()
                // Every GEX row of a library shares species and chemistry; main checks samplesheets agree.
                .map { lid, metas, runs ->
                    def m = metas[0]
                    [lid, [library_id: lid, run_name: m.run_name], velocity_chemistry(m.chemistry),
                     file(m.species == 'human' ? params.spliceu_index_human : params.spliceu_index_mouse),
                     runs.toSorted { r -> r[0] }.collectMany { r -> r[1].toSorted { f -> f.name } }]
                }
                .join(ch_cellbender_bc.map { meta, bc -> [meta.library_id, bc] })
                .map { row -> row.tail() }

            SIMPLEAF_VELOCITY(ch_velocity_input)
        }
}
