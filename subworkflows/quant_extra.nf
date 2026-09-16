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
                    def meta   = metas[0] + [library_id: library_id]
                    def bam    = file("${outs}/unassigned_alignments.bam")
                    def bai    = file("${outs}/unassigned_alignments.bam.bai")
                    def wl     = file(get_viral_whitelist(meta.chemistry, params.tenx_barcodes_dir))
                    def schem  = get_simpleaf_chemistry(meta.chemistry)
                    [meta, bam, bai, wl, schem]
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
                .map { meta, fastq_dir, _fqs ->
                    [ meta.library_id, meta, fastq_dir, get_velocity_chemistry(meta.chemistry) ]
                }
                .groupTuple(by: 0)
                .map { library_id, metas, fastq_dirs, chems ->
                    def sorted_metas = metas.toSorted { m -> m.toString() }
                    def uniq_chems   = chems.toUnique()
                    assert uniq_chems.size() == 1 :
                        "library ${library_id} has conflicting velocity chemistries: ${uniq_chems}"
                    def meta = sorted_metas[0] + [library_id: library_id, run_name: run_name]
                    [ library_id, meta, fastq_dirs.toUnique().toSorted().join(','), uniq_chems[0] ]
                }
                .join(
                    ch_cellbender_bc.map { meta, bc -> [ meta.library_id, bc ] },
                    by: 0
                )
                .multiMap { _library_id, meta, fastq_dirs, chemistry, barcodes ->
                    def is_human = meta.species == 'human'
                    def idx = file(is_human ? params.spliceu_index_human : params.spliceu_index_mouse)
                    input: [ meta, fastq_dirs, chemistry, barcodes ]
                    index: idx
                }
                .set { ch_velocity_split }

            SIMPLEAF_VELOCITY(
                ch_velocity_split.input,
                ch_velocity_split.index
            )
        }
}
