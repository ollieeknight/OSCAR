include { MULTI_CONFIG }      from '../modules/count_gex'
include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas_list, fastq_dirs, adt_csv]
    main:
        ch_libraries
            .map { lid, metas, dirs, adt_csv ->
                def ml       = (metas as List).findAll { it != null }.sort { it.id }
                def meta     = ml.find { it.modality == 'GEX' } ?: ml[0]
                def is_human = meta.species == 'human'
                def ref_gex  = is_human ? params.ref_human : params.ref_mouse
                def ref_vdj  = is_human ? params.ref_vdj_human : params.ref_vdj_mouse
                def has_vdj  = ml.any { it.modality in ['VDJ-T', 'VDJ-B'] }
                def has_adt  = ml.any { it.modality in ['ADT', 'HTO'] }
                def create_bam = (is_human && meta.n_donors > 1) ? 'true' : 'false'

                def lines = ['[gene-expression]',
                             "reference,${ref_gex}",
                             "create-bam,${create_bam}"]
                if (meta.assay in ['DOGMA', 'Multiome'])      lines << 'chemistry,ARC-v1'
                else if (meta.assay == 'Flex' && meta.chemistry) lines << "chemistry,${meta.chemistry}"

                if (has_vdj)  lines += ['', '[vdj]', "reference,${ref_vdj}"]
                if (has_adt && adt_csv.name != 'NO_FILE')
                    lines += ['', '[feature]', "reference,${adt_csv.toAbsolutePath()}"]

                lines += ['', '[libraries]', 'fastq_id,fastqs,feature_types']

                // Per (dir × modality): count reads in the R1 file to exclude BCL Convert
                // placeholder FASTQs. Placeholders can exceed 10 MB but have <100 real reads.
                // Threshold matches cellranger's own auto-detection minimum (10 000 reads) —
                // any flowcell contributing fewer reads than this would cause cellranger to
                // fail with TXRNGR10001 during chemistry auto-detection.
                def min_reads = 10000
                def lib_checks = (dirs as List).sort().collectMany { dir ->
                    ml.collect { m ->
                        def ft = m.modality == 'GEX'          ? 'Gene Expression'      :
                                 m.modality in ['ADT','HTO']  ? 'Antibody Capture'     :
                                 m.modality == 'VDJ-T'        ? 'VDJ-T'                :
                                 m.modality == 'VDJ-B'        ? 'VDJ-B'                :
                                 m.modality == 'CRISPR'       ? 'CRISPR Guide Capture' : 'Gene Expression'
                        [
                            "r1=\$(find \"${dir}\" -maxdepth 2 -name \"${m.id}*_R1_*.fastq.gz\" -print -quit 2>/dev/null)",
                            "n_reads=0; [ -n \"\$r1\" ] && n_reads=\$(zcat \"\$r1\" 2>/dev/null | head -n ${min_reads * 4} | awk 'NR%4==1' | wc -l)",
                            "[ \"\$n_reads\" -ge ${min_reads} ] && echo \"${m.id},${dir},${ft}\" >> multi_config.csv || true"
                        ].join('\n')
                    }
                }

                [lid, ml, lines.join('\n'), lib_checks.join('\n'), adt_csv]
            }
            .set { ch_ready }

        MULTI_CONFIG(ch_ready)
        CELLRANGER_MULTI(MULTI_CONFIG.out.config)
    emit:
        outs     = CELLRANGER_MULTI.out.outs
        versions = CELLRANGER_MULTI.out.versions
}