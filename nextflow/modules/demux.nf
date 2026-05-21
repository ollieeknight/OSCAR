// ─── Base mask / OverrideCycles helpers ───────────────────────────────────────
// Ported from functions.sh:100-321 (check_base_masks_step1/2/3)
// get_base_mask returns [cellranger_cmd, filter_option, mask_string]
// get_override_cycles converts mask_string → BCL Convert OverrideCycles format

def get_base_mask(assay, chemistry, index_type, modality, num_reads) {
    def masks3 = [
        'SI_SC3Pv2_GEX':        'Y26n*,I8n*,Y98n*',
        'SI_SC3Pv2_ADT':        'Y26n*,I8n*,Y98n*',
        'SI_SC3Pv2_HTO':        'Y26n*,I8n*,Y98n*',
        'SI_SC3Pv3_GEX':        'Y28n*,I8n*,Y90n*',
        'SI_SC3Pv3_ADT':        'Y28n*,I8n*,Y90n*',
        'SI_SC3Pv3_HTO':        'Y28n*,I8n*,Y90n*',
        'SI_SC5P_GEX':          'Y26n*,I8n*,Y90n*',
        'SI_SC5P_ADT':          'Y26n*,I8n*,Y90n*',
        'SI_SC5P_HTO':          'Y26n*,I8n*,Y90n*',
        'SI_SC5P_VDJ':          'Y26n*,I8n*,Y90n*',
        'SI_DOGMA_ARCv1_ADT':   'Y24n*,I8n*,Y90n*',
        'SI_DOGMA_ARCv1_HTO':   'Y28n*,I8n*,Y90n*',
    ]
    def masks4 = [
        'SI_SC3Pv2_GEX':        'Y26n*,I8n*,N*,Y98n*',
        'SI_SC3Pv2_ADT':        'Y26n*,I8n*,N*,Y98n*',
        'SI_SC3Pv2_HTO':        'Y26n*,I8n*,N*,Y98n*',
        'SI_SC3Pv3_GEX':        'Y28n*,I8n*,N*,Y90n*',
        'SI_SC3Pv3_ADT':        'Y28n*,I8n*,N*,Y90n*',
        'SI_SC3Pv3_HTO':        'Y28n*,I8n*,N*,Y90n*',
        'SI_SC5P_GEX':          'Y26n*,I8n*,N*,Y90n*',
        'SI_SC5P_ADT':          'Y26n*,I8n*,N*,Y90n*',
        'SI_SC5P_HTO':          'Y26n*,I8n*,N*,Y90n*',
        'SI_SC5P_VDJ':          'Y26n*,I8n*,N*,Y90n*',
        'SI_DOGMA_ARCv1_ADT':   'Y24n*,I8n*,N*,Y90n*',
        'SI_DOGMA_ARCv1_HTO':   'Y28n*,I8n*,N*,Y90n*',
        'DI_SC3Pv2_GEX':        'Y26n*,I8n*,N*,Y98n*',
        'DI_SC3Pv2_ADT':        'Y26n*,I8n*,N*,Y98n*',
        'DI_SC3Pv2_HTO':        'Y26n*,I8n*,N*,Y98n*',
        'DI_SC3Pv3_GEX':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC3Pv3_ADT':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC3Pv3_HTO':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC3Pv4_GEX':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC3Pv4_ADT':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC3Pv4_HTO':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC5P_GEX':          'Y26n*,I10n*,I10n*,Y90n*',
        'DI_SC5P_ADT':          'Y26n*,I10n*,I10n*,Y90n*',
        'DI_SC5P_HTO':          'Y26n*,I10n*,I10n*,Y90n*',
        'DI_SC5P_VDJ':          'Y26n*,I10n*,I10n*,Y90n*',
        'DI_SC5Pv3_GEX':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC5Pv3_ADT':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC5Pv3_HTO':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_SC5Pv3_VDJ':        'Y28n*,I10n*,I10n*,Y90n*',
        'DI_Multiome_ARCv1_GEX':  'Y28n*,I10n*,I10n*,Y90n*',
        'DI_Multiome_ARCv1_ATAC': '50n*,I8n*,Y24n*,Y49n*',
        'DI_DOGMA_ARCv1_GEX':   'Y28n*,I10n*,I10n*,Y90n*',
        'DI_DOGMA_ARCv1_ATAC':  'Y100n*,I8n*,Y24n*,Y100n*',
        'DI_DOGMA_ARCv1_ADT':   'Y28n*,I8n*,N*,Y90n*',
        'DI_DOGMA_ARCv1_HTO':   'Y28n*,I8n*,N*,Y90n*',
        'DI_ATAC_ATAC':         'Y50n*,I8n*,Y16n*,Y50n*',
        'DI_ASAP_ATAC':         'Y100n*,I8n*,Y16n*,Y100n*',
        'DI_ASAP_ADT':          'Y100n*,I8n*,Y16n*,Y100n*',
        'DI_ASAP_HTO':          'Y100n*,I8n*,Y16n*,Y100n*',
        'DI_ASAP_GENO':         'Y100n*,I8n*,Y16n*,Y100n*',
    ]

    def masks = (num_reads == 3) ? masks3 : masks4
    def chem  = chemistry.replaceAll(/[-_ ]/, '')

    def key
    if (assay in ['CITE', 'GEX']) {
        if      (chem.startsWith('SC3Pv2'))                     key = "${index_type}_SC3Pv2_${modality}"
        else if (chem.startsWith('SC3Pv3'))                     key = "${index_type}_SC3Pv3_${modality}"
        else if (chem.startsWith('SC3Pv4'))                     key = "${index_type}_SC3Pv4_${modality}"
        else if (chem.startsWith('SC5P') && chem.contains('v3')) key = "${index_type}_SC5Pv3_${modality}"
        else if (chem.startsWith('SC5P'))                       key = "${index_type}_SC5P_${modality}"
        else key = null
    } else if (assay == 'Multiome') {
        key = "DI_Multiome_ARCv1_${modality}"
    } else if (assay == 'DOGMA') {
        key = (modality == 'ATAC') ? "DI_DOGMA_ARCv1_ATAC" : "${index_type}_DOGMA_ARCv1_${modality}"
    } else if (assay == 'ASAP') {
        key = "DI_ASAP_${modality}"
    } else if (assay == 'ATAC') {
        key = "DI_ATAC_ATAC"
    } else {
        key = null
    }

    def base_mask = key ? masks.get(key) : null
    if (!base_mask)
        error "Cannot determine base mask: assay=${assay} chem=${chemistry} index_type=${index_type} modality=${modality} num_reads=${num_reads} (key=${key})"
    return base_mask
}

// Convert bcl2fastq mask format → BCL Convert OverrideCycles format.
// Example: 'Y28n*,I10n*,I10n*,Y90n*' → 'Y28N*;I10N*;I10N*;Y90N*'
// index_seqs: resolved [is_dual, rows] from main.nf — used to fix SI-on-DI runs
// where the base mask has I10N* for both index positions but the library is
// single-index (8bp). On 4-read flow cells the I1 position is clamped to the
// actual sequence length and the I2 position is fully masked (N*).
// ATAC cell-barcode reads in position 3 are Y reads, not I, so they are untouched.
def get_override_cycles(assay, chemistry, index_type, modality, num_reads, index_seqs = null) {
    def mask
    try {
        mask = get_base_mask(assay, chemistry, index_type, modality, num_reads)
    } catch (e) {
        return null  // no entry for this num_reads; caller falls back to the other read count
    }

    def oc = mask.replace(',', ';').replace('n', 'N')

    if (num_reads == 4 && index_seqs != null && !index_seqs.is_dual) {
        def parts     = oc.split(';') as List
        def seq_len   = index_seqs.rows[0].i7.length()
        def i_indices = (0..<parts.size()).findAll { parts[it].startsWith('I') }
        def fixed     = (0..<parts.size()).collect { idx ->
            if      (idx == i_indices[0])                               "I${seq_len}N*"
            else if (i_indices.size() > 1 && idx == i_indices[1])      'N*'
            else                                                         parts[idx]
        }
        return fixed.join(';')
    }

    return oc
}

// ─── GENERATE_SAMPLESHEET ────────────────────────────────────────────────────
// Builds a BCL Convert V2 SampleSheet for one demux group.
// Reads RunInfo.xml at runtime to get actual cycle lengths, then resolves
// OverrideCycles wildcards (*) to exact counts required by BCL Convert 4.x.
// Input channel: [demux_key, metas_list, bcl_dir]

process GENERATE_SAMPLESHEET {
    tag "$demux_key"
    label 'process_low'
    container "${params.container_bclconvert}"

    input:
    tuple val(demux_key), val(metas), path(bcl_dir)

    output:
    tuple val(demux_key), val(metas), path(bcl_dir), path("SampleSheet.csv"), emit: samplesheet

    script:
    def meta        = metas[0]
    def oc_4        = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 4, meta.index_seqs)
    def oc_3        = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 3, meta.index_seqs) ?: oc_4
    def is_dual     = metas.any { m -> m.index_seqs.is_dual }
    def data_header = is_dual ? 'Sample_ID,Index,Index2' : 'Sample_ID,Index'
    def data_rows   = metas.collectMany { m ->
        m.index_seqs.rows.collect { row ->
            is_dual ? "${m.id},${row.i7},${row.get('i5', '')}" : "${m.id},${row.i7}"
        }
    }.join('\n')

    """
    mapfile -t cycles < <(grep -o 'NumCycles="[0-9]*"' ${bcl_dir}/RunInfo.xml | grep -o '[0-9]*')
    num_reads=\${#cycles[@]}
    r1=\${cycles[0]}
    i1=\${cycles[1]}
    if [ "\$num_reads" -eq 4 ]; then
        i2=\${cycles[2]}
        r2=\${cycles[3]}
        raw_oc="${oc_4}"
        read_lens=("\$r1" "\$i1" "\$i2" "\$r2")
    else
        r2=\${cycles[2]}
        raw_oc="${oc_3}"
        read_lens=("\$r1" "\$i1" "\$r2")
    fi

    # Resolve * to exact counts — BCL Convert 4.x rejects wildcards
    IFS=';' read -ra oc_parts <<< "\$raw_oc"
    expanded=()
    for i in "\${!oc_parts[@]}"; do
        part="\${oc_parts[\$i]}"
        len="\${read_lens[\$i]}"
        if [[ "\$part" == *'*' ]]; then
            base="\${part%\\*}"
            used=\$(echo "\$base" | grep -oE '[0-9]+' | awk '{s+=\$1}END{print s+0}')
            rest=\$((len - used))
            if [ "\$rest" -gt 0 ]; then
                expanded+=("\${base}\${rest}")
            else
                expanded+=("\$(echo "\$base" | sed 's/[A-Z]\$//')")
            fi
        else
            expanded+=("\$part")
        fi
    done
    override_cycles=\$(IFS=';'; echo "\${expanded[*]}")

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
        echo "OverrideCycles,\$override_cycles"
        echo 'BarcodeMismatchesIndex1,1'
        [ "${is_dual}" = "true" ] && echo 'BarcodeMismatchesIndex2,1' || true
        echo ''
        echo '[BCLConvert_Data]'
    } > SampleSheet.csv

    cat >> SampleSheet.csv << 'DATAEOF'
${data_header}
${data_rows}
DATAEOF
    """
}

// ─── BCL_TO_FASTQ ─────────────────────────────────────────────────────────────
// Runs BCL Convert using a pre-built SampleSheet from GENERATE_SAMPLESHEET.
// Input channel: [demux_key, metas_list, bcl_dir, samplesheet]

process BCL_TO_FASTQ {
    tag "$demux_key"
    label 'process_medium'   // overridden to 16c/32GB/12h via withName: 'BCL_TO_FASTQ'
    container "${params.container_bclconvert}"
    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${params.outdir}/${run}_fastq"
    }, mode: 'copy', pattern: 'fastqs/*.fastq.gz',
        saveAs: { fn -> fn.tokenize('/')[-1] }
    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${params.outdir}/${run}_fastq"
    }, mode: 'copy', pattern: 'fastqs/Reports/Top_Unknown_Barcodes*.csv',
        saveAs: { fn -> fn.tokenize('/')[-1] }

    input:
    tuple val(demux_key), val(metas), path(bcl_dir), path(samplesheet)

    output:
    tuple val(metas), val(bcl_dir.name), path("fastqs/*.fastq.gz"), emit: fastqs
    path "fastqs/Reports/Top_Unknown_Barcodes*.csv",                 optional: true, emit: unknown_barcodes
    path "versions.yml",                                              emit: versions

    script:
    def has_gex = metas.any { it.modality == 'GEX' }
    def num_unknown = has_gex ? '50' : '0'
    """
    bcl-convert \\
        --bcl-input-directory              ${bcl_dir} \\
        --output-directory                 fastqs \\
        --sample-sheet                     ${samplesheet} \\
        --no-lane-splitting                true \\
        --force \\
        --bcl-only-matched-reads           true \\
        --bcl-num-parallel-tiles           4 \\
        --bcl-num-conversion-threads       ${task.cpus} \\
        --bcl-enable-tile-metrics          false \\
        --bcl-enable-adapter-cycle-metrics false \\
        --num-unknown-barcodes-reported    ${num_unknown}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl-convert: \$(bcl-convert --version 2>&1 | head -1 | sed 's/BCL Convert //')
    END_VERSIONS
    """
}

// ─── FALCO ────────────────────────────────────────────────────────────────────
// Per-file QC. One job per R-read FASTQ (R1/R2/R3); index reads (I1/I2) skipped.

process FALCO {
    tag "$fastq_name"
    label 'process_low'
    container "${params.container_falco}"
    publishDir { "${params.outdir}/${run_name}_fastq/falco" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_name), path(fastq)

    output:
    path "${fastq_name}/", emit: report
    path "versions.yml",   emit: versions

    script:
    """
    mkdir -p ${fastq_name}
    falco -t ${task.cpus} ${fastq} -o ${fastq_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: \$(falco --version 2>&1 | sed 's/falco //')
    END_VERSIONS
    """
}

// ─── MULTIQC ──────────────────────────────────────────────────────────────────

process MULTIQC {
    label 'process_low'
    container "${params.container_multiqc}"
    publishDir { "${params.outdir}/${params.run_name}_outs/multiqc" }, mode: 'copy'

    input:
    path(reports)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data
    path "versions.yml",        emit: versions

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    multiqc ${config} --force -o . .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
