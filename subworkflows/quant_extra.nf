include { VIRAL_DETECT; SIMPLEAF_VELOCITY } from '../modules/quant_extra'
include { get_viral_whitelist; get_simpleaf_chemistry; get_velocity_chemistry } from '../lib/chemistry'


workflow QUANT_EXTRA {
    take:
        ch_gex_outs
        ch_gex_fastqs
        ch_cellbender_bc
        extras
        run_name

    main:
        if ('viral' in extras) {
            ch_gex_outs
                .filter { _library_id, metas, _outs -> metas[0].species == 'human' }
                .map { library_id, metas, outs ->
                    def chem = metas[0].chemistry
                    [[library_id: library_id, run_name: metas[0].run_name],
                     file("${outs}/unassigned_alignments.bam"), file("${outs}/unassigned_alignments.bam.bai"),
                     file(get_viral_whitelist(chem, params.tenx_barcodes_dir)), get_simpleaf_chemistry(chem)]
                }
                .set { ch_viral_input }

            VIRAL_DETECT(
                ch_viral_input,
                file(params.viral_piscem_index),
                file(params.viral_t2g),
                file(params.bamtofastq_bin)
            )
        }

        if ('velocity' in extras) {
            ch_gex_fastqs
                .filter { meta, _fastq_dir, _fqs ->
                    meta.modality == 'GEX' && get_velocity_chemistry(meta.chemistry) != null
                }
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, [fastq_dir, fqs.findAll { f -> f.name =~ /_R[12]_/ }]] }
                .groupTuple(by: 0)
                // Every GEX row of a library shares species and chemistry; main checks samplesheets agree.
                .map { library_id, metas, runs ->
                    def m = metas[0]
                    [library_id, [library_id: library_id, run_name: run_name], get_velocity_chemistry(m.chemistry),
                     file(m.species == 'human' ? params.spliceu_index_human : params.spliceu_index_mouse),
                     runs.toSorted { r -> r[0] }.collectMany { r -> r[1].toSorted { f -> f.name } }]
                }
                .join(ch_cellbender_bc.map { meta, bc -> [meta.library_id, bc] })
                .map { row -> row.tail() }
                .set { ch_velocity_input }

            SIMPLEAF_VELOCITY(ch_velocity_input)
        }
}
