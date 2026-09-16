include { CELLRANGER_MQC } from '../modules/demux'
include { MULTIQC }        from '../modules/fastq_qc'


workflow REPORT {
    take:
        ch_gex_outs
        ch_fastp_reports
        ch_demux_mqc

    main:
        ch_gex_outs
            .map { _library_id, metas, outs ->
                [metas.toSorted { m -> m.id }[0].run_name, outs]
            }
            .groupTuple(by: 0)
            .map { run_name, outs -> [run_name, outs.toSorted { o -> o.toString() }] }
            .set { ch_cellranger_outs }

        CELLRANGER_MQC(ch_cellranger_outs)

        ch_fastp_reports
            .join(ch_demux_mqc, by: [0, 1], remainder: true)
            .map { run_name, fastq_dir, reports, demux_mqc ->
                [run_name, fastq_dir, (reports ?: []) + (demux_mqc ? [demux_mqc] : [])]
            }
            .map { run_name, fastq_dir, files -> [run_name, [fastq_dir, files]] }
            .join(CELLRANGER_MQC.out.mqc, by: 0, remainder: true)
            .map { run_name, pair, cr_mqc ->
                def fastq_dir = pair ? pair[0] : null
                def files     = pair ? pair[1] : []
                [run_name, fastq_dir, files + (cr_mqc ? [cr_mqc] : [])]
            }
            .filter { _run_name, fastq_dir, files -> fastq_dir && files }
            .set { ch_multiqc_in }

        MULTIQC(ch_multiqc_in)
}
