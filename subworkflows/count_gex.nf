include { MULTI_CONFIG }      from '../modules/count_gex'
include { CELLRANGER_MULTI } from '../modules/count_gex'

workflow COUNT_GEX {
    take:
        ch_libraries   // [library_id, metas_list, fastq_dirs, adt_csv]
    main:
        ch_libraries
            .map { lid, metas, dirs, adt_csv ->
                def ml       = (metas as List).findAll { it != null }
                def meta     = ml[0]
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

                // Build library check shell lines
                def lib_checks = (dirs as List).collectMany { dir ->
                    ml.collect { m ->
                        def ft = m.modality == 'GEX'          ? 'Gene Expression'      :
                                 m.modality in ['ADT','HTO']  ? 'Antibody Capture'     :
                                 m.modality == 'VDJ-T'        ? 'VDJ-T'                :
                                 m.modality == 'VDJ-B'        ? 'VDJ-B'                :
                                 m.modality == 'CRISPR'       ? 'CRISPR Guide Capture' : 'Gene Expression'
                        "find \"${dir}\" -maxdepth 2 -name \"${m.id}*.fastq.gz\" " +
                        "-print -quit 2>/dev/null | grep -q . && " +
                        "echo \"${m.id},${dir},${ft}\" >> multi_config.csv || true"
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