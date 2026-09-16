include { get_chemistry_family } from '../lib/chemistry'

def get_override_cycles(assay, chemistry, index_type, modality, num_reads, index_seqs = null, index_len = 10) {
    def masks = [
        'SI_SC3Pv2_GEX':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv2_ADT':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv2_HTO':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv3_GEX':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv3_ADT':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv3_HTO':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4_GEX':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4_ADT':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4_HTO':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC5P_GEX':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_ADT':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_HTO':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_VDJ':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_ADT':   [3: 'Y24N*;I8N*;Y90N*',           4: 'Y24N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_HTO':   [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_SC3Pv2_GEX':        [4: 'Y26N*;I8N*;N*;Y98N*'],
        'DI_SC3Pv2_ADT':        [4: 'Y26N*;I8N*;N*;Y98N*'],
        'DI_SC3Pv2_HTO':        [4: 'Y26N*;I8N*;N*;Y98N*'],
        'DI_SC3Pv3_GEX':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv3_ADT':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv3_HTO':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv4_GEX':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv4_ADT':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv4_HTO':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P_GEX':          [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P_ADT':          [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P_HTO':          [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P_VDJ':          [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3_GEX':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3_ADT':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3_HTO':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3_VDJ':        [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_Multiome_ARCv1_GEX':  [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_Multiome_ARCv1_ATAC': [4: '50N*;I8N*;Y24N*;Y49N*'],
        'DI_DOGMA_ARCv1_GEX':   [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_DOGMA_ARCv1_ATAC':  [4: 'Y100N*;I8N*;Y24N*;Y100N*'],
        'DI_DOGMA_ARCv1_ADT':   [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_DOGMA_ARCv1_HTO':   [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_ATAC_ATAC':         [4: 'Y50N*;I8N*;Y16N*;Y50N*'],
        'DI_ASAP_ATAC':         [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_ADT':          [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_HTO':          [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_Flex-v2_GEX':       [4: 'Y*;I10;I10;Y*'],
    ]

    def mod_key = modality.replaceAll(/^VDJ-[TB]$/, 'VDJ').replaceAll(/^CRISPR$/, 'GEX')
    def key
    if (assay in ['CITE', 'GEX']) {
        def family = get_chemistry_family(chemistry)
        key = (family in ['SC3Pv2', 'SC3Pv3', 'SC3Pv4', 'SC5P', 'SC5Pv3'])
            ? "${index_type}_${family}_${mod_key}"
            : null
    } else if (assay == 'Flex') {
        key = "DI_Flex-v2_GEX"
    } else if (assay == 'Multiome') {
        key = "DI_Multiome_ARCv1_${mod_key}"
    } else if (assay == 'DOGMA') {
        key = (modality == 'ATAC') ? "DI_DOGMA_ARCv1_ATAC" : "${index_type}_DOGMA_ARCv1_${mod_key}"
    } else if (assay == 'ASAP') {
        key = "DI_ASAP_${mod_key}"
    } else if (assay == 'ATAC') {
        key = "DI_ATAC_ATAC"
    } else {
        key = null
    }

    if (!key || !masks.containsKey(key))
        error "Cannot determine OverrideCycles: assay=${assay} chem=${chemistry} index_type=${index_type} modality=${modality} num_reads=${num_reads} (key=${key})"

    def mask_entry = masks[key]
    if (!mask_entry.containsKey(num_reads))
        return null

    def oc = mask_entry[num_reads]

    if (index_len == 8) {
        oc = oc.replaceAll(/I10N\*/, 'I8N2')
    }

    if (num_reads == 4 && index_seqs != null && !index_seqs.is_dual) {
        oc = apply_si_on_di_correction(oc, index_seqs.rows[0].i7.length())
    }

    return oc
}

def apply_si_on_di_correction(String oc, Integer seq_len, Integer index1_cycles = 10) {
    def parts     = oc.split(';') as List
    def i_indices = (0..<parts.size()).findAll { i -> parts[i].startsWith('I') }

    if (i_indices.isEmpty()) return oc

    def fixed = (0..<parts.size()).collect { idx ->
        if (idx == i_indices[0]) {
            def remaining = index1_cycles - seq_len
            remaining > 0 ? "I${seq_len}N${remaining}" : "I${seq_len}"
        } else if (i_indices.size() > 1 && idx == i_indices[1]) {
            'N*'
        } else {
            parts[idx]
        }
    }

    return fixed.join(';')
}

process GENERATE_SAMPLESHEET {
    tag "$demux_key"
    container "${params.container_bclconvert}"

    input:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), val(is_dual), val(sample_specs)

    output:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), path("SampleSheet.csv"), emit: samplesheet

    script:
    def specs = sample_specs.collect { sp ->
        def oc4 = get_override_cycles(sp.assay, sp.chemistry, sp.index_type, sp.modality, 4, [is_dual: sp.is_dual, rows: [[i7: sp.i7]]], sp.index_len)
        def oc3 = get_override_cycles(sp.assay, sp.chemistry, sp.index_type, sp.modality, 3, [is_dual: sp.is_dual, rows: [[i7: sp.i7]]], sp.index_len) ?: oc4
        [id: sp.id, i7: sp.i7, i5: sp.i5 ?: '', is_dual: sp.is_dual, oc4: oc4, oc3: oc3]
    }
    def per_sample = specs.collect { sp -> "${sp.oc4}|${sp.oc3}" }.unique().size() > 1

    def spec_tsv = specs.collect { sp -> "${sp.id}\t${sp.i7}\t${sp.i5 ?: '-'}\t${sp.oc4}\t${sp.oc3}" }.join('\n')

    """
    mapfile -t cycles < <(grep -o 'NumCycles="[0-9]*"' ${bcl_dir}/RunInfo.xml | grep -o '[0-9]*')
    num_reads=\${#cycles[@]}
    r1=\${cycles[0]}
    i1=\${cycles[1]}
    if [ "\$num_reads" -eq 4 ]; then
        i2=\${cycles[2]}
        r2=\${cycles[3]}
        read_lens=("\$r1" "\$i1" "\$i2" "\$r2")
    else
        r2=\${cycles[2]}
        read_lens=("\$r1" "\$i1" "\$r2")
    fi

    cat > sample_specs.tsv << 'SPECEOF'
${spec_tsv}
SPECEOF

    # Per-row expanded masks, in samplesheet row order.
    : > row_masks.txt
    while IFS=\$'\\t' read -r sid si7 si5 soc4 soc3; do
        [ -z "\$sid" ] && continue
        if [ "\$num_reads" -eq 4 ]; then raw="\$soc4"; else raw="\$soc3"; fi
        expand_override_cycles.sh "\$raw" "\${read_lens[@]}" >> row_masks.txt || exit 1
    done < sample_specs.tsv

    # Global OverrideCycles only when every row agrees; per-sample column
    # otherwise. Never both — bcl-convert rejects a setting given twice.
    per_sample=${per_sample}
    override_cycles=\$(head -n1 row_masks.txt)

    {
        echo '[Header]'
        echo 'FileFormatVersion,2'
        echo ''
        echo '[Reads]'
        echo "Read1Cycles,\$r1"
        echo "Index1Cycles,\$i1"
        [ "\$num_reads" -eq 4 ] && echo "Index2Cycles,\$i2" || true
        echo "Read2Cycles,\$r2"
        echo ''
        echo '[BCLConvert_Settings]'
        [ "\$per_sample" = "true" ] || echo "OverrideCycles,\$override_cycles"
        echo 'BarcodeMismatchesIndex1,1'
        if [ "${is_dual}" = "true" ] && [ "\$per_sample" != "true" ]; then
            echo 'BarcodeMismatchesIndex2,1'
        fi
        echo ''
        echo '[BCLConvert_Data]'
    } > SampleSheet.csv

    {
        header='Sample_ID,Index'
        [ "${is_dual}" = "true" ] && header="\$header,Index2"
        [ "\$per_sample" = "true" ] && header="\$header,OverrideCycles"
        echo "\$header"
        paste sample_specs.tsv row_masks.txt | while IFS=\$'\\t' read -r sid si7 si5 soc4 soc3 mask; do
            [ -z "\$sid" ] && continue
            [ "\$si5" = "-" ] && si5=""
            row="\$sid,\$si7"
            [ "${is_dual}" = "true" ] && row="\$row,\$si5"
            [ "\$per_sample" = "true" ] && row="\$row,\$mask"
            echo "\$row"
        done
    } >> SampleSheet.csv
    """
}

process CLEAN_FASTQ_DIR {
    tag "$bcl_name"
    executor 'local'

    input:
    tuple val(bcl_name), val(bcl_parent)

    output:
    tuple val(bcl_name), val(true), emit: done

    script:
    def run = bcl_name.replaceAll(/_bcl.*$/, '')
    """
    mkdir -p "${bcl_parent}/${run}_fastq"
    rm -f "${bcl_parent}/${run}_fastq"/*.fastq.gz
    """
}


process BCLCONVERT {
    tag "$demux_key L$lane"
    container "${params.container_bclconvert}"
    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${bcl_parent}/${run}_fastq"
    }, mode: 'copy', pattern: "fastqs/*.fastq.gz", saveAs: { fn -> file(fn).name }

    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${bcl_parent}/${run}_fastq/Reports/${demux_key}_L${lane}"
    }, mode: 'copy', pattern: "fastqs/Reports/*", saveAs: { fn -> file(fn).name }

    input:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), path(samplesheet), val(lane)

    output:
    tuple val(demux_key), val(metas), val(bcl_dir.name), val(bcl_parent), path("fastqs/*.fastq.gz"), emit: fastqs
    tuple val(demux_key), val(metas), val(bcl_dir.name), val(bcl_parent), val(lane), path("fastqs/Reports/Demultiplex_Stats.csv"), path("fastqs/Reports/Top_Unknown_Barcodes.csv"), emit: reports

    script:
    def n_tiles      = Math.max(1, (task.cpus / 8).toInteger())
    def n_convert    = Math.max(1, (task.cpus / 8).toInteger())
    def n_compress   = Math.max(1, (task.cpus / 2).toInteger())
    def n_decompress = Math.max(1, (task.cpus / 4).toInteger())
    """
    rm -rf fastqs/

    bcl-convert \\
        --bcl-input-directory              ${bcl_dir} \\
        --output-directory                 fastqs \\
        --sample-sheet                     ${samplesheet} \\
        --bcl-only-lane                    ${lane} \\
        --bcl-num-parallel-tiles           ${n_tiles} \\
        --bcl-num-conversion-threads       ${n_convert} \\
        --bcl-num-compression-threads      ${n_compress} \\
        --bcl-num-decompression-threads    ${n_decompress}
    """
}

process DEMUX_QC {
    tag "$run_name"
    container "${params.container_multiqc}"
    publishDir { "${fastq_dir}" }, mode: 'copy', pattern: "*.csv"

    input:
    tuple val(run_name), val(fastq_dir), path(stats, stageAs: 'stats/*/'), path(unknown, stageAs: 'unknown/*/')
    path indexes_dir

    output:
    tuple val(run_name), val(fastq_dir), path("${run_name}_demux_summary.csv"),  emit: summary
    tuple val(run_name), val(fastq_dir), path("${run_name}_demux_warnings.csv"), emit: warnings
    tuple val(run_name), val(fastq_dir), path("${run_name}_demux_mqc.csv"), emit: mqc

    script:
    """
    # Concatenate the per-group/per-lane Reports, keeping a single header.
    awk 'FNR==1 && NR!=1 { next } { print }' stats/*/*.csv   > all_demultiplex_stats.csv
    awk 'FNR==1 && NR!=1 { next } { print }' unknown/*/*.csv > all_top_unknown.csv

    demux_qc.py \\
        --stats    all_demultiplex_stats.csv \\
        --unknown  all_top_unknown.csv \\
        --indexes  ${indexes_dir} \\
        --prefix   ${run_name} \\
        --dropout-ratio ${params.demux_qc_dropout_ratio} \\
        --unknown-pct   ${params.demux_qc_unknown_pct}

    # MultiQC custom-content: header comment block makes it a named section.
    {
        echo '# id: oscar_demux'
        echo '# section_name: Demultiplexing (bcl-convert)'
        echo '# description: "Reads per library from bcl-convert Demultiplex_Stats.csv, as a percent of each lane total (Undetermined included)."'
        echo '# plot_type: table'
        cat ${run_name}_demux_summary.csv
    } > ${run_name}_demux_mqc.csv
    """
}

process CELLRANGER_MQC {
    tag "$run_name"
    container "${params.container_multiqc}"

    input:
    tuple val(run_name), path(outs, stageAs: 'outs/*/')

    output:
    tuple val(run_name), path("${run_name}_cellranger_mqc.csv"), emit: mqc, optional: true

    script:
    """
    cellranger_mqc.py --root outs --out ${run_name}_cellranger_mqc.csv
    """
}
