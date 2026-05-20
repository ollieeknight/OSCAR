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

// ─── BCL_TO_FASTQ ─────────────────────────────────────────────────────────────
// Direct BCL Convert (bclconvert) demultiplexing.
// One job per unique {assay}_{index_type}_{chemistry}_{modality} combination.
// Generates V2 SampleSheet from resolved meta.index_seqs at runtime.
// Input channel: [demux_key, metas_list, bcl_dir]

process BCL_TO_FASTQ {
    tag "$demux_key"
    label 'process_medium'   // overridden to 16c/32GB/12h via withName: 'BCL_TO_FASTQ'
    container "${params.container_bclconvert}"

    input:
    tuple val(demux_key), val(metas), path(bcl_dir)

    output:
    tuple val(metas), path("fastqs/*.fastq.gz"), emit: fastqs
    path "versions.yml",                          emit: versions

    script:
    def meta     = metas[0]
    def oc_4     = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 4, meta.index_seqs)
    def oc_3     = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 3, meta.index_seqs) ?: oc_4
    def is_dual  = metas.any { m -> m.index_seqs.is_dual }

    // Build [BCLConvert_Data] section (static — known from samplesheet)
    def data_header = is_dual ? 'Sample_ID,Index,Index2' : 'Sample_ID,Index'
    def data_rows   = metas.collectMany { m ->
        m.index_seqs.rows.collect { row ->
            is_dual ? "${m.id},${row.i7},${row.get('i5', '')}" : "${m.id},${row.i7}"
        }
    }.join('\n')

    """
    # Parse RunInfo.xml: count reads and extract cycle lengths (pure bash, no python)
    mapfile -t cycles < <(grep -o 'NumCycles="[0-9]*"' ${bcl_dir}/RunInfo.xml | grep -o '[0-9]*')
    num_reads=\${#cycles[@]}
    r1=\${cycles[0]}
    i1=\${cycles[1]}
    if [ "\$num_reads" -eq 4 ]; then
        i2=\${cycles[2]}
        r2=\${cycles[3]}
        override_cycles="${oc_4}"
    else
        r2=\${cycles[2]}
        override_cycles="${oc_3}"
    fi

    # Write V2 SampleSheet header + settings (dynamic, depends on RunInfo.xml)
    {
        echo '[Header]'
        echo 'FileFormatVersion,2'
        echo ''
        echo '[Reads]'
        echo "Read1Cycles,\$r1"
        echo "Index1Cycles,\$i1"
        [ "\$num_reads" -eq 4 ] && echo "Index2Cycles,\$i2" || true
        echo "Read2Cycles,\$([ "\$num_reads" -eq 4 ] && echo \$r2 || echo \$r2)"
        echo ''
        echo '[BCLConvert_Settings]'
        echo "OverrideCycles,\$override_cycles"
        echo 'BarcodeMismatchesIndex1,1'
        [ "${is_dual}" = "true" ] && echo 'BarcodeMismatchesIndex2,1' || true
        echo ''
        echo '[BCLConvert_Data]'
    } > SampleSheet.csv

    # Append sample rows (injected by Groovy, not shell variables)
    cat >> SampleSheet.csv << 'DATAEOF'
${data_header}
${data_rows}
DATAEOF

    # Run BCL Convert
    bcl-convert \\
        --bcl-input-directory   ${bcl_dir} \\
        --output-directory      fastqs \\
        --sample-sheet          SampleSheet.csv \\
        --bcl-num-parallel-tiles ${task.cpus} \\
        --no-lane-splitting     false

    # Remove undetermined reads
    rm -f fastqs/Undetermined_*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl-convert: \$(bcl-convert --version 2>&1 | head -1 | sed 's/BCL Convert //')
    END_VERSIONS
    """
}

// ─── FALCO ────────────────────────────────────────────────────────────────────
// Batch QC on all FASTQs from one demux job.

process FALCO {
    tag "$demux_key"
    label 'process_medium'
    container "${params.container_falco}"

    input:
    tuple val(demux_key), path(fastqs)

    output:
    path "falco_${demux_key}/", emit: report
    path "versions.yml",        emit: versions

    script:
    """
    mkdir -p falco_${demux_key}
    find . -maxdepth 1 -name '*.fastq.gz' \\
        | xargs -P ${task.cpus} -I{} falco -t 1 {} -o falco_${demux_key}

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
    publishDir "${params.outdir}/multiqc", mode: 'copy'

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
