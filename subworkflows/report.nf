include { MULTIQC } from '../modules/demux'

workflow REPORT {
    take:
        ch_reports   // [run_name, [report_dirs]]

    main:
        MULTIQC(ch_reports)

    emit:
        report   = MULTIQC.out.report
}
