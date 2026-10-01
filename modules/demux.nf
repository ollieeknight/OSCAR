include { get_chemistry_family } from '../lib/chemistry'

// CITE/GEX masks depend only on index type and chemistry family; VDJ needs a 5' chemistry.
def get_override_cycles(Map sp, int num_reads) {
    def masks = [
        'SI_SC3Pv2':              [3: 'Y26N*;I8N*;Y98N*',  4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv3':              [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4':              [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC5P':                [3: 'Y26N*;I8N*;Y90N*',  4: 'Y26N*;I8N*;N*;Y90N*'],
        'DI_SC3Pv2':              [4: 'Y26N*;I8N*;N*;Y98N*'],
        'DI_SC3Pv3':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC3Pv4':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_SC5P':                [4: 'Y26N*;I10N*;I10N*;Y90N*'],
        'DI_SC5Pv3':              [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'SI_DOGMA_ARCv1_ADT':     [3: 'Y24N*;I8N*;Y90N*',  4: 'Y24N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_HTO':     [3: 'Y28N*;I8N*;Y90N*',  4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_Multiome_ARCv1_GEX':  [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_Multiome_ARCv1_ATAC': [4: '50N*;I8N*;Y24N*;Y49N*'],
        'DI_DOGMA_ARCv1_GEX':     [4: 'Y28N*;I10N*;I10N*;Y90N*'],
        'DI_DOGMA_ARCv1_ATAC':    [4: 'Y100N*;I8N*;Y24N*;Y100N*'],
        'DI_DOGMA_ARCv1_ADT':     [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_DOGMA_ARCv1_HTO':     [4: 'Y28N*;I8N*;N*;Y90N*'],
        'DI_ATAC_ATAC':           [4: 'Y50N*;I8N*;Y16N*;Y50N*'],
        'DI_ASAP_ATAC':           [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_ADT':            [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_ASAP_HTO':            [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
        'DI_Flex-v2_GEX':         [4: 'Y*;I10;I10;Y*'],
    ]

    def mod = sp.modality.replaceAll(/^VDJ-[TB]$/, 'VDJ').replaceAll(/^CRISPR$/, 'GEX')
    def family = sp.assay in ['CITE', 'GEX'] ? get_chemistry_family(sp.chemistry) : null
    def key = family && (mod in ['GEX', 'ADT', 'HTO'] || (mod == 'VDJ' && family.startsWith('SC5P'))) ? "${sp.index_type}_${family}"
            : sp.assay in ['CITE', 'GEX'] ? null
            : sp.assay == 'Flex'          ? 'DI_Flex-v2_GEX'
            : sp.assay == 'Multiome'      ? "DI_Multiome_ARCv1_${mod}"
            : sp.assay == 'DOGMA'         ? (sp.modality == 'ATAC' ? 'DI_DOGMA_ARCv1_ATAC' : "${sp.index_type}_DOGMA_ARCv1_${mod}")
            : sp.assay == 'ASAP'          ? "DI_ASAP_${mod}"
            : sp.assay == 'ATAC'          ? 'DI_ATAC_ATAC'
            : null
    def oc = masks[key?.toString()]
    if (!oc)
        error "Cannot determine OverrideCycles: assay=${sp.assay} chem=${sp.chemistry} index_type=${sp.index_type} modality=${sp.modality} (key=${key})"
    oc = oc[num_reads]
    if (!oc) return null

    if (sp.index_len == 8) oc = oc.replaceAll(/I10N\*/, 'I8N2')
    if (num_reads == 4 && !sp.is_dual) oc = apply_si_on_di_correction(oc, sp.index_len)
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
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), val(is_dual), val(sample_specs)

    output:
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), path("SampleSheet.csv"), emit: samplesheet

    script:
    def specs = sample_specs.collect { sp ->
        def oc4 = get_override_cycles(sp, 4)
        def oc3 = get_override_cycles(sp, 3) ?: oc4
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
    # otherwise. Never both: bcl-convert rejects a setting given twice.
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

// Clears FASTQs from an earlier demux so the published folder holds only this run's files.
process CLEAN_FASTQ_DIR {
    tag "$fastq_dir"
    executor 'local'

    input:
    val(fastq_dir)

    output:
    val(fastq_dir), emit: done

    script:
    """
    mkdir -p "${fastq_dir}"
    rm -f "${fastq_dir}"/*.fastq.gz
    """
}

process BCLCONVERT {
    tag "$demux_key L$lane"
    container "${params.container_bclconvert}"
    publishDir { fastq_dir }, mode: 'copy', pattern: "fastqs/*.fastq.gz", saveAs: { fn -> file(fn).name }
    publishDir { "${fastq_dir}/Reports/${demux_key}_L${lane}" }, mode: 'copy', pattern: "fastqs/Reports/*", saveAs: { fn -> file(fn).name }

    input:
    tuple val(demux_key), path(bcl_dir), val(fastq_dir), path(samplesheet), val(lane)

    output:
    tuple val(demux_key), val(fastq_dir), path("fastqs/*.fastq.gz"), emit: fastqs
    tuple val(fastq_dir), path("fastqs/Reports/Demultiplex_Stats.csv"), path("fastqs/Reports/Top_Unknown_Barcodes.csv"), emit: reports

    script:
    def n_tiles      = Math.max(1, (task.cpus / 8).toInteger())
    def n_compress   = Math.max(1, (task.cpus / 2).toInteger())
    def n_decompress = Math.max(1, (task.cpus / 4).toInteger())
    """
    bcl-convert \\
        --bcl-input-directory              ${bcl_dir} \\
        --output-directory                 fastqs \\
        --sample-sheet                     ${samplesheet} \\
        --bcl-only-lane                    ${lane} \\
        --bcl-num-parallel-tiles           ${n_tiles} \\
        --bcl-num-conversion-threads       ${n_tiles} \\
        --bcl-num-compression-threads      ${n_compress} \\
        --bcl-num-decompression-threads    ${n_decompress}
    """
}

process DEMUX_QC {
    tag "$run_name"
    container "${params.container_python}"
    publishDir { "${fastq_dir}" }, mode: 'copy', pattern: "*.csv"

    input:
    tuple val(run_name), val(fastq_dir), path(stats, stageAs: 'stats/*/'), path(unknown, stageAs: 'unknown/*/')
    path indexes_dir

    output:
    tuple val(run_name), val(fastq_dir), path("${run_name}_demux_summary.csv"),  emit: summary
    tuple val(run_name), val(fastq_dir), path("${run_name}_demux_warnings.csv"), emit: warnings
    tuple val(run_name), val(fastq_dir), path("${run_name}_flow_cell_overview.csv"), emit: flowcell_overview

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
    """
}
