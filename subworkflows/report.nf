include { CELLRANGER_MQC } from '../modules/demux'
include { MULTIQC }        from '../modules/fastq_qc'

// Run-level reporting. Fans in from every counting branch, so it runs after
// counting: one report covers demultiplexing (fastp + the demux summary table)
// and Cell Ranger metrics.

workflow REPORT {
    take:
        ch_gex_outs        // [library_id, metas, outs/] from COUNT_GEX
        ch_fastp_reports   // [run_name, fastq_dir, [report_files]] from DEMUX
        ch_demux_mqc       // [run_name, fastq_dir, demux_mqc.csv] from DEMUX

    main:
        // cellranger multi writes per_sample_outs/*/metrics_summary.csv per
        // library; MultiQC's built-in cellranger module cannot read a `multi`
        // web summary, so CELLRANGER_MQC pivots the CSVs into a table.
        ch_gex_outs
            .map { _library_id, metas, outs -> [metas[0].run_name, outs] }
            .groupTuple(by: 0)
            .set { ch_cellranger_outs }

        CELLRANGER_MQC(ch_cellranger_outs)

        // remainder:true on both joins: a run may finish demux with no
        // countable library, and MultiQC should still report the demux.
        ch_fastp_reports
            .join(ch_demux_mqc, by: [0, 1], remainder: true)
            .map { run_name, fastq_dir, reports, demux_mqc ->
                [run_name, fastq_dir, (reports ?: []) + (demux_mqc ? [demux_mqc] : [])]
            }
            .map { run_name, fastq_dir, files -> [run_name, [fastq_dir, files]] }
            .join(CELLRANGER_MQC.out.mqc, by: 0, remainder: true)
            .map { run_name, pair, cr_mqc ->
                // remainder:true pads the missing side with null.
                def fastq_dir = pair ? pair[0] : null
                def files     = pair ? pair[1] : []
                [run_name, fastq_dir, files + (cr_mqc ? [cr_mqc] : [])]
            }
            .filter { _run_name, fastq_dir, files -> fastq_dir && files }
            .set { ch_multiqc_in }

        MULTIQC(ch_multiqc_in)
}
