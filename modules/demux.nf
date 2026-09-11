include { get_chemistry_family } from '../lib/chemistry'

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
        'SI_SC3Pv4_GEX':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4_ADT':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
        'SI_SC3Pv4_HTO':        [3: 'Y28N*;I8N*;Y90N*',           4: 'Y28N*;I8N*;N*;Y90N*'],
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
        // GEM-X Flex v2 (Fixed RNA Profiling): full cycles — N-masking causes bcl-convert cbcl stall
        'DI_Flex-v2_GEX':       [4: 'Y*;I10;I10;Y*'],
    ]

    // Resolve key from assay/chemistry/index_type/modality.
    // Chemistry family comes from the registry in lib/chemistry.nf so that
    // adding a chemistry there is enough to make it resolvable here.
    // Normalise modality: VDJ-T/VDJ-B → VDJ (same read structure), CRISPR → GEX (same read structure)
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
// index1_cycles is the flow cell's actual Index1 cycle count from RunInfo.xml.
// It defaults to 10 only because that is what this was hardcoded to before the
// per-row path existed; pass the real value whenever it is known, or an 8bp
// index on a 12-cycle index read produces a mask 2 cycles short of the read.
def apply_si_on_di_correction(String oc, Integer seq_len, Integer index1_cycles = 10) {
    def parts     = oc.split(';') as List
    def i_indices = (0..<parts.size()).findAll { parts[it].startsWith('I') }

    if (i_indices.isEmpty()) return oc  // no index positions, no correction needed

    def fixed = (0..<parts.size()).collect { idx ->
        if (idx == i_indices[0]) {
            // First index: use actual seq_len (typically 8bp), mask remainder
            def remaining = index1_cycles - seq_len  // e.g., 10 - 8 = 2
            remaining > 0 ? "I${seq_len}N${remaining}" : "I${seq_len}"
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
//
// Mixed-index groups (8bp TruSeq SI alongside 10bp 10x DI in one lane) are
// demultiplexed in a SINGLE bcl-convert call using a per-sample OverrideCycles
// column in [BCLConvert_Data], supported since BCL Convert 4.1.5 (in use here:
// 4.5.4). Per Illumina's documentation a setting may be specified globally OR
// per-sample but NOT both, so when the per-sample column is emitted the
// OverrideCycles line is omitted from [BCLConvert_Settings].
//
// This replaces an earlier attempt that put 8bp Index values under a single
// global 10bp mask. That does not work: bcl-convert matches the Index column
// against the cycle count declared by OverrideCycles exactly and does not
// prefix-match, so it aborts with
//   "has an index of length 8 bases, but a length of 10 was expected".
// The 8bp rows instead get their own mask that uses 8 index cycles and masks
// the remaining ones (I8N2), with i5 fully masked (N10) for single-index rows.
//
// Input channel: [demux_key, metas_list, bcl_dir, bcl_parent, is_dual, sample_specs]
// sample_specs: one map per samplesheet row —
//   [id, i7, i5, index_len, is_dual, assay, chemistry, index_type, modality]

process GENERATE_SAMPLESHEET {
    tag "$demux_key"
    container "${params.container_bclconvert}"

    input:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), val(is_dual), val(sample_specs)

    output:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), path("SampleSheet.csv"), emit: samplesheet

    script:
    // Resolve each row's mask template up front (Groovy side), keeping the
    // runtime wildcard expansion in the shell where RunInfo.xml is readable.
    // A group is "mixed" when its rows do not all share one mask template and
    // index length — that is exactly when the per-sample column is required.
    def specs = sample_specs.collect { sp ->
        def oc4 = get_override_cycles(sp.assay, sp.chemistry, sp.index_type, sp.modality, 4, [is_dual: sp.is_dual, rows: [[i7: sp.i7]]], sp.index_len)
        def oc3 = get_override_cycles(sp.assay, sp.chemistry, sp.index_type, sp.modality, 3, [is_dual: sp.is_dual, rows: [[i7: sp.i7]]], sp.index_len) ?: oc4
        [id: sp.id, i7: sp.i7, i5: sp.i5 ?: '', is_dual: sp.is_dual, oc4: oc4, oc3: oc3]
    }
    def per_sample = specs.collect { "${it.oc4}|${it.oc3}" }.unique().size() > 1

    // TSV consumed by the shell loop below: id, i7, i5, 4-read mask, 3-read mask.
    // '-' stands in for an absent i5: bash `read` collapses empty fields even
    // with IFS=$'\t', which shifts every subsequent column. The placeholder is
    // turned back into an empty string when the row is written.
    def spec_tsv = specs.collect { "${it.id}\t${it.i7}\t${it.i5 ?: '-'}\t${it.oc4}\t${it.oc3}" }.join('\n')

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

    # Resolve * to exact counts — BCL Convert 4.x rejects wildcards.
    # Emits the expanded mask on stdout; fails loudly if a mask needs more
    # cycles than the read actually has.
    expand_mask() {
        local raw="\$1"
        local -a oc_parts expanded
        IFS=';' read -ra oc_parts <<< "\$raw"
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
                elif [ "\$rest" -eq 0 ]; then
                    expanded+=("\${base%[A-Z]}")
                else
                    echo "ERROR: read \$((i+1)) has only \${len} cycles but mask '\${part}' requires at least \${used} cycles. Check RunInfo.xml and OverrideCycles mask." >&2
                    return 1
                fi
            else
                expanded+=("\$part")
            fi
        done
        (IFS=';'; echo "\${expanded[*]}")
    }

    cat > sample_specs.tsv << 'SPECEOF'
${spec_tsv}
SPECEOF

    # Per-row expanded masks, in samplesheet row order.
    : > row_masks.txt
    while IFS=\$'\\t' read -r sid si7 si5 soc4 soc3; do
        [ -z "\$sid" ] && continue
        if [ "\$num_reads" -eq 4 ]; then raw="\$soc4"; else raw="\$soc3"; fi
        expand_mask "\$raw" >> row_masks.txt || exit 1
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
        # BarcodeMismatchesIndex2 is only valid when index2 is used for
        # demultiplexing across ALL samples. On a mixed sheet the single-index
        # rows mask i5 entirely, so the setting is omitted there.
        if [ "${is_dual}" = "true" ] && [ "\$per_sample" != "true" ]; then
            echo 'BarcodeMismatchesIndex2,1'
        fi
        echo ''
        echo '[BCLConvert_Data]'
    } > SampleSheet.csv

    # Data section. Index2 column present whenever any row in the group is
    # dual-indexed; single-index rows leave it blank and mask i5 via their own
    # OverrideCycles entry.
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

// ─── CLEAN_FASTQ_DIR ──────────────────────────────────────────────────────────
// Wipes previously-published *.fastq.gz for one bcl_dir before any of this
// run's BCLCONVERT tasks publish into the same directory. Runs once per unique
// bcl_dir, gated ahead of BCLCONVERT via combine(by:0) in the subworkflow.
// publishDir only ever adds/overwrites files — it never removes a file that a
// prior invocation published under a since-changed sample sheet (different
// bcl-convert-assigned S-number), so without this, old and new S-numbered
// fastq for the same library silently coexist after a re-run.

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

// ─── BCLCONVERT ─────────────────────────────────────────────────────────────

process BCLCONVERT {
    tag "${demux_key}_L${lane}"
    container "${params.container_bclconvert}"
    // FASTQs are published beside their source flowcell (bcl_parent), not under
    // params.outdir — so each run's reads land in its own directory when several
    // flowcells are merged via --extra_bcl_dirs.
    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${bcl_parent}/${run}_fastq"
    }, mode: 'copy', pattern: "fastqs/*.fastq.gz", saveAs: { fn -> file(fn).name }
    // bcl-convert writes Demultiplex_Stats.csv / Quality_Metrics.csv /
    // Top_Unknown_Barcodes.csv into <output-directory>/Reports for free. Keep
    // them: without these, diagnosing a bad demux means decompressing the
    // multi-GB Undetermined FASTQ by hand. One Reports dir per demux group and
    // lane, so key the destination by both to avoid groups overwriting
    // each other's stats.
    publishDir {
        def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
        "${bcl_parent}/${run}_fastq/Reports/${demux_key}_L${lane}"
    }, mode: 'copy', pattern: "fastqs/Reports/*", saveAs: { fn -> file(fn).name }

    input:
    tuple val(demux_key), val(metas), path(bcl_dir), val(bcl_parent), path(samplesheet), val(lane)

    output:
    tuple val(demux_key), val(metas), val(bcl_dir.name), path("fastqs/*.fastq.gz"), emit: fastqs
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

// ─── DEMUX_QC ─────────────────────────────────────────────────────────────────
// Summarises bcl-convert's own demultiplexing stats for one run and flags a
// suspect demux. Each BCLCONVERT task (one per demux group per lane) writes its
// own Reports/, so all of a run's Demultiplex_Stats.csv are concatenated here
// before analysis -- otherwise each group is judged only against itself and a
// whole missing group looks normal.
//
// Advisory only: warnings never fail the run, since the FASTQ are still valid.

process DEMUX_QC {
    tag "$run_name"
    container "${params.container_multiqc}"
    publishDir { "${fastq_dir}" }, mode: 'copy', pattern: "*.csv"

    input:
    // Every demux group writes a Reports/ dir with identically-named files, so
    // they are staged into numbered subdirectories ('stats/1/', 'stats/2/'...).
    // A flat 'stats/*' collides the moment a run has more than one group.
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

// ─── CELLRANGER_MQC ───────────────────────────────────────────────────────────
// Pivots `cellranger multi` per-sample metrics into a MultiQC custom-content
// table. MultiQC's built-in cellranger module matches web summaries on
// '"subcommand":"count"' / '"subcommand":"vdj"' only, so a `cellranger multi`
// web_summary.html is never picked up and the metrics go unreported.

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

// ─── FASTP ────────────────────────────────────────────────────────────────────
// Per-file QC, report only. One job per R-read FASTQ (R1/R2/R3); index reads
// (I1/I2) skipped. No -o/--out1, so fastp writes reports and no filtered reads.

process FASTP {
    tag "$fastq_name"
    container "${params.container_fastp}"
    // Published beside the source flowcell's FASTQs (fastq_dir), not under
    // params.outdir — otherwise every --extra_bcl_dirs run drops a stray
    // {run}_fastq/fastp tree into the primary run's directory.
    publishDir { "${fastq_dir}/fastp" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), val(fastq_name), path(fastq)

    output:
    tuple val(run_name), val(fastq_dir), path("${run_name}_${fastq_name}.{json,html}"), emit: report
    // Completion signal keyed by the FASTQ this task verified. FASTQ_QC uses it
    // to release the file downstream, so a fastp failure blocks its library.
    tuple val(fastq_name), val(true), emit: checked

    script:
    """
    fastp \\
        --in1                       ${fastq} \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --thread                    ${task.cpus} \\
        --json                      ${run_name}_${fastq_name}.json \\
        --html                      ${run_name}_${fastq_name}.html
    """
}

// ─── MULTIQC ──────────────────────────────────────────────────────────────────

process MULTIQC {
    container "${params.container_multiqc}"
    // Beside the source flowcell's FASTQs, matching FASTP.
    publishDir { "${fastq_dir}/multiqc" }, mode: 'copy'

    input:
    tuple val(run_name), val(fastq_dir), path(reports)

    output:
    path "multiqc_report.html",      emit: report
    path "multiqc_report_data/",     emit: data

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    multiqc ${config} --force --filename multiqc_report -o . .
    """
}

// ─── VALIDATE_FASTQ ───────────────────────────────────────────────────────────
// Lightweight validation step that runs gzip -t on each individual fastq file.
// Fully distributed across Slurm nodes and benefits from Nextflow caching.

process VALIDATE_FASTQ {
    tag "$meta.id"
    container "${params.container_pigz}"

    input:
    tuple val(meta), val(fastq_dir), path(fastq), val(fastq_name)

    output:
    tuple val(meta), val(fastq_dir), path(fastq), emit: fastq

    script:
    """
    pigz -t -f -p ${task.cpus} ${fastq}
    """
}
