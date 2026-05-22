// ─── Override Cycles Helpers ──────────────────────────────────────────────────
// BCL Convert OverrideCycles for each assay/chemistry/index/modality combination.
// Format: Y=data read, I=index, N=masked cycle. Semicolons separate read segments.
// Example: 'Y28N*;I10N*;I10N*;Y90N*' = 28bp read1, 10bp index1, 10bp index2, 90bp read2.
// SI-on-DI (single-index on 4-read dual-index flow cell): detected at runtime and corrected.

def get_override_cycles(assay, chemistry, index_type, modality, num_reads, index_seqs = null, index_len = 10) {
    // Masks in BCL Convert format (native). 3-read flow cells use position 2; 4-read add index2.
    // Stored for 10bp indices (10x kits). For 8bp indices (TruSeq), I10 → I8N2 substitution applied below.
    def masks = [
        // 3-read (SI only, position 2 is R2)
        'SI_SC3Pv2_GEX':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv2_ADT':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv2_HTO':        [3: 'Y26N*;I8N*;Y98N*',           4: 'Y26N*;I8N*;N*;Y98N*'],
        'SI_SC3Pv3_GEX':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv3_ADT':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv3_HTO':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC5P_GEX':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_ADT':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_HTO':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_SC5P_VDJ':          [3: 'Y26N*;I8N*;Y90N*',           4: 'Y26N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_ADT':   [3: 'Y24N*;I8N*;Y90N*',           4: 'Y24N*;I8N*;N*;Y90N*'],
        'SI_DOGMA_ARCv1_HTO':   [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        // 4-read dual-index (DI only, position 2/3 are I2/R2)
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
        'DI_ASAP_GENO':         [4: 'Y100N*;I8N*;Y16N*;Y100N*'],
    ]

    // Resolve key from assay/chemistry/index_type/modality
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

    if (!key || !masks.containsKey(key))
        error "Cannot determine OverrideCycles: assay=${assay} chem=${chemistry} index_type=${index_type} modality=${modality} num_reads=${num_reads} (key=${key})"

    def mask_entry = masks[key]
    if (!mask_entry.containsKey(num_reads))
        return null  // no entry for this num_reads; caller falls back to the other read count

    def oc = mask_entry[num_reads]

    // Adjust OverrideCycles for 8bp indices (TruSeq): replace I10 with I8N2, keep I8N* as-is
    if (index_len == 8) {
        oc = oc.replaceAll(/I10N\*/, 'I8N2')
    }

    // SI-on-DI correction: single-index library on 4-read dual-index flow cell.
    // Clamp index1 to actual sequence length, fully mask index2.
    if (num_reads == 4 && index_seqs != null && !index_seqs.is_dual) {
        oc = apply_si_on_di_correction(oc, index_seqs.rows[0].i7.length())
    }

    return oc
}

// Correct OverrideCycles for single-index (8bp) library on 4-read dual-index flow cell.
// Input: 'Y28N*;I10N*;I10N*;Y90N*', seq_len=8
// Output: 'Y28N*;I8N2;N*;Y90N*' (I1 uses 8bp+2masked, I2 fully masked)
def apply_si_on_di_correction(String oc, Integer seq_len) {
    def parts     = oc.split(';') as List
    def i_indices = (0..<parts.size()).findAll { parts[it].startsWith('I') }

    if (i_indices.isEmpty()) return oc  // no index positions, no correction needed

    def fixed = (0..<parts.size()).collect { idx ->
        if (idx == i_indices[0]) {
            // First index: use actual seq_len (typically 8bp), mask remainder
            def remaining = 10 - seq_len  // e.g., 10 - 8 = 2
            "I${seq_len}N${remaining > 0 ? remaining : ''}"
        } else if (i_indices.size() > 1 && idx == i_indices[1]) {
            // Second index: fully masked (no i5 present)
            'N*'
        } else {
            // Data reads (Y) and cell barcodes (ATAC): untouched
            parts[idx]
        }
    }

    return fixed.join(';')
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
    def index_len   = meta.index_seqs?.rows[0]?.i7?.length() ?: 10
    // Verify all metas in group have same index length (should be guaranteed by grouping key, but be defensive)
    if (metas.any { m -> (m.index_seqs?.rows[0]?.i7?.length() ?: 10) != index_len }) {
        error "Demux group has mixed index lengths: ${metas.collect { m -> m.id + '=' + (m.index_seqs?.rows[0]?.i7?.length() ?: 10) }.join(', ')}"
    }
    def oc_4        = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 4, meta.index_seqs, index_len)
    def oc_3        = get_override_cycles(meta.assay, meta.chemistry, meta.index_type, meta.modality, 3, meta.index_seqs, index_len) ?: oc_4
    def is_dual     = metas.any { m -> m.index_seqs.is_dual }
    def data_header = is_dual ? 'Sample_ID,Index,Index2' : 'Sample_ID,Index'
    def data_rows   = metas.collectMany { m ->
        m.index_seqs.rows.collect { row ->
            is_dual ? "${m.id},${row.i7},${row.get('i5', '')}" : "${m.id},${row.i7}"        }
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
    }, mode: 'copy', pattern: 'fastqs/Reports/Top_Unknown_Barcodes_*.csv',
        saveAs: { fn -> fn.tokenize('/')[-1] }
    input:
    tuple val(demux_key), val(metas), path(bcl_dir), path(samplesheet)

    output:
    tuple val(metas), val(bcl_dir.name), path("fastqs/*.fastq.gz"), emit: fastqs
    path "fastqs/Reports/Top_Unknown_Barcodes_*.csv",                optional: true, emit: unknown_barcodes
    path "versions.yml",                                              emit: versions

    script:
    def modality    = metas[0].modality
    def num_unknown = modality == 'GEX' ? '50' : '0'
    """
    bcl-convert \\
        --bcl-input-directory              ${bcl_dir} \\
        --output-directory                 fastqs \\
        --sample-sheet                     ${samplesheet} \\
        --no-lane-splitting                true \\
        --force \\
        --bcl-num-parallel-tiles           4 \\
        --bcl-num-conversion-threads       ${task.cpus} \\
        --bcl-enable-tile-metrics          false \\
        --bcl-enable-adapter-cycle-metrics false \\
        --bcl-only-matched-reads           true \\
        --num-unknown-barcodes-reported    ${num_unknown}

    for f in fastqs/Reports/Top_Unknown_Barcodes.csv fastqs/Reports/Top_Unknown_Barcodes_L*.csv; do
        [ -f "\$f" ] && mv "\$f" "\${f/Top_Unknown_Barcodes/Top_Unknown_Barcodes_${modality}}" || true
    done

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
