# PLAN: Viral Detection in OSCAR

Map single-cell FASTQ files to human-infecting viral sequences to detect viral transcripts per library or per cell.

---

## Context

- Input: short-read Illumina FASTQs from 10x GEX/ATAC/CITE/DOGMA assays
- Current pipeline routes everything to cellranger → QC; no viral awareness
- Goal: detect and quantify human-infecting viral reads, ideally with single-cell resolution for GEX

---

## Phase 0: Documentation Discovery

Before implementing any approach, gather:

1. **Tool availability** — check biocontainers tags for Kraken2, Viral-Track, STARsolo viral reference builds
2. **Database sources** — confirm download URLs and update cadence for RVDB, NCBI RefSeq viral, ViPR
3. **Existing nf-core modules** — check nf-core/modules for `kraken2/kraken2`, `star/align`, and viral ref-building modules to avoid reinventing
4. **Viral-Track paper/repo** — confirm STARsolo compatibility and output format (https://github.com/PieroMastroberardino/Viral-Track)

Allowed references:
- RVDB nucleotide: https://rvdb.dbi.udel.edu/ (curated human-infecting, ~100k sequences)
- NCBI RefSeq viral complete genomes: accessible via `datasets download virus` CLI
- nf-core/modules repo for Kraken2/STAR module signatures

---

## Method Comparison

| Method | Resolution | Speed | Known-virus only | SC-native | Integration effort |
|---|---|---|---|---|---|
| A. Custom cellranger reference (viral + human) | per-cell UMI | slowest | yes | yes (native) | low |
| B. Viral-Track (STAR + viral db) | per-cell barcode | medium | semi (large db) | yes | medium |
| C. Kraken2 on unmapped reads | per-library (reads) | fastest | no (k-mer) | no | low |
| D. STAR alignment to viral ref (unmapped reads only) | per-library | medium | yes | no | medium |
| E. De novo assembly (Trinity/SPAdes) of unmapped reads | per-library | slow | no | no | high |

**Recommendation**: implement A + C as complementary approaches:
- A for targeted known-virus per-cell quantification (EBV, CMV, HHV-6, HTLV-1, HIV, etc.)
- C for unbiased screening of all reads using Kraken2 + RVDB/RefSeq

---

## Method A: Custom Cellranger Reference (Known Viruses)

### Strategy
Concatenate curated viral genomes + GTF annotations to the existing human reference before calling `cellranger mkref`. Cellranger then counts viral UMIs per cell alongside human genes automatically.

### Steps

**A1 — Build viral reference set**
- Select target viruses: EBV (HHV-4), CMV (HHV-5), HHV-6A/B, HHV-7, HHV-8/KSHV, HIV-1/2, HTLV-1/2, SARS-CoV-2, Influenza A/B, RSV, Adenovirus, HPV (high-risk: 16/18/31/33/45)
- Download: NCBI RefSeq via `datasets download virus genome taxon <taxid> --complete-only --filename viral.zip`
- Alternatively: download from ViPR (https://www.viprbrc.org) for curated sets
- Format: extract FASTA + GFF3 → convert GFF3 to GTF using `gffread`

**A2 — Mask low-complexity / human-homologous viral regions**
- Run `dustmasker` on viral FASTA to soft-mask repetitive regions
- Optional: BLAST viral seqs against human genome, mask regions with >80% identity over >50 bp
- Prevents false-positive viral reads from human repetitive elements

**A3 — Concatenate to human reference**
- Prepend `virus_` prefix to all viral chromosome names to avoid contig name collision
- Concatenate `hg38.fa + viral.fa` → `hg38_viral.fa`
- Concatenate `hg38.gtf + viral.gtf` → `hg38_viral.gtf`
- Run `cellranger mkref --genome=hg38_viral --fasta=hg38_viral.fa --genes=hg38_viral.gtf`
- This is a one-time reference build step, stored at `params.ref_human_viral`

**A4 — Pipeline integration**
- Add `params.ref_human_viral = null` (default: off)
- In COUNT_GEX subworkflow: if `params.ref_human_viral` set, pass it to cellranger multi instead of `params.ref_human`
- No new processes needed — viral counts appear in `filtered_feature_bc_matrix` under `Gene Expression` feature type
- Viral gene names prefixed `virus_*` are trivially filterable downstream

**A5 — Verification**
```
grep "^virus_" {outs}/filtered_feature_bc_matrix/features.tsv.gz | wc -l   # should > 0
```

### Pros / Cons
- Pro: zero extra compute, cell-level resolution, uses existing cellranger infrastructure
- Con: reference rebuild required per viral set update; only detects pre-selected viruses; inflates reference size slightly

---

## Method B: Viral-Track (STARsolo per-cell)

### Strategy
Viral-Track aligns raw FASTQs to a large viral genome database using STAR (with cell barcode/UMI awareness), producing per-cell viral read counts. Designed specifically for 10x Chromium GEX.

### Steps

**B1 — Build Viral-Track STAR index**
- Download RVDB-nt (`C-RVDBvXX.fasta.gz`, human-virus curated): shttps://rvdb.dbi.udel.edu/
- Build STAR index: `STAR --runMode genomeGenerate --genomeDir viral_star_index --genomeFastaFiles rvdb.fa --genomeSAindexNbases 7`  (small genome → use `--genomeSAindexNbases 7`)
- Store at `params.viral_star_index`

**B2 — New process: VIRAL_TRACK**
```nextflow
process VIRAL_TRACK {
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
    label 'process_high'
    input:
        tuple val(meta), val(fastq_dirs), path(viral_index)
    output:
        tuple val(meta), path("${meta.library_id}_viral/"), emit: counts
        path "versions.yml",                                 emit: versions
    script:
    """
    STAR --soloType CB_UMI_Simple \
         --soloCBwhitelist ${params.gex_whitelist} \
         --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
         --readFilesIn <(find ${fastq_dirs.join(' ')} -name '*_R2_*.fastq.gz' | sort | tr '\\n' ',') \
                       <(find ${fastq_dirs.join(' ')} -name '*_R1_*.fastq.gz' | sort | tr '\\n' ',') \
         --readFilesCommand zcat \
         --genomeDir ${viral_index} \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --soloFeatures Gene \
         --outFileNamePrefix ${meta.library_id}_viral/ \
         --runThreadN ${task.cpus}
    ...
    """
}
```

**B3 — Pipeline integration**
- Add `params.viral_star_index = null` (default: off)
- In main.nf: after DEMUX/from_fastq, fork `ch_routed.gex` into VIRAL_TRACK when param set
- Viral counts land at `{library_id}_viral/Solo.out/Gene/`

**B4 — Verification**
- Check that viral barcode matrix has non-zero entries for positive controls

### Pros / Cons
- Pro: broad viral discovery (RVDB covers ~99% human-infecting viruses); cell-level resolution
- Con: RVDB alignment produces many spurious hits; requires careful post-filtering; STAR index large (~20 GB)

---

## Method C: Kraken2 Taxonomic Classification (Screening)

### Strategy
Classify all reads (or unmapped reads from cellranger BAM) using Kraken2 + a human-infecting viral database. Fast, per-library, no alignment required.

### Steps

**C1 — Build Kraken2 viral database**
Option 1 (standard): `kraken2-build --download-library viral --db kraken2_viral_db && kraken2-build --build --db kraken2_viral_db`
Option 2 (RVDB-based): download RVDB FASTA → `kraken2-build --add-to-library rvdb.fa --db kraken2_rvdb_db && kraken2-build --build --db kraken2_rvdb_db`
Option 3 (human-infecting only): use `kraken2-build` with `--taxid` filter list of human-infecting viral taxids from NCBI

Recommended: Option 1 (standard NCBI viral RefSeq) for speed; Option 2 for sensitivity. Store at `params.kraken2_viral_db`.

**C2 — Input: raw FASTQs or unmapped BAM**
- Raw FASTQs (simpler): classify all reads → human reads dominate, but Kraken2 handles this with `--confidence 0.1`
- Unmapped reads (lower FP): `samtools view -f 4 cellranger_bam.bam | samtools fastq | kraken2 ...`

For OSCAR: start with raw FASTQs (no dependency on cellranger completing first).

**C3 — New process: KRAKEN2_VIRAL**
```nextflow
process KRAKEN2_VIRAL {
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'
    label 'process_high'
    input:
        tuple val(meta), path(fastqs)
        path db
    output:
        tuple val(meta), path("${meta.library_id}_kraken2_report.txt"), emit: report
        path "versions.yml",                                              emit: versions
    script:
    """
    kraken2 --db ${db} \
            --threads ${task.cpus} \
            --report ${meta.library_id}_kraken2_report.txt \
            --gzip-compressed \
            --confidence 0.1 \
            --paired ${fastqs.findAll { it.name =~ /_R1_/ }.join(' ')} \
                     ${fastqs.findAll { it.name =~ /_R2_/ }.join(' ')}
    ...
    """
}
```

**C4 — Output/reporting**
- Kraken2 report is compatible with MultiQC (`multiqc --module kraken`)
- Add Kraken2 reports to `MULTIQC` call in main.nf
- Bracken can re-estimate abundance at species level from the report if needed

**C5 — Verification**
- Spike in known viral reads (e.g. SARS-CoV-2 synthetic reads) and confirm classification
- Check MultiQC kraken section renders

### Pros / Cons
- Pro: fastest option (minutes per library); no index rebuild for new viruses; integrates with MultiQC
- Con: per-library only (no cell resolution); k-mer approach has lower specificity; reports reads not UMIs

---

## Method D: STAR Alignment to Viral Reference (Unmapped Reads Only)

### Strategy
Extract unmapped reads from cellranger BAM, then align to a viral STAR index. More specific than Kraken2 (alignment-based), less noisy than aligning all reads.

### Steps

**D1 — Extract unmapped reads**
```bash
samtools view -f 4 -b possorted_genome_bam.bam | samtools fastq -1 R1.fq.gz -2 R2.fq.gz -
```

**D2 — Align to viral STAR index** (same index as Method B)
```bash
STAR --genomeDir viral_star_index \
     --readFilesIn R2.fq.gz R1.fq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --runThreadN 16
```

**D3 — Pipeline integration**
- Runs after COUNT_GEX completes (depends on `{outs_dir}/possorted_genome_bam.bam`)
- Output: per-library viral read counts by gene (ReadsPerGene.out.tab)

### Pros / Cons
- Pro: cleanest signal (human reads already removed); alignment-based specificity
- Con: depends on cellranger completing first (longer total wall-time); no cell-level resolution

---

## Method E: De Novo Assembly of Unmapped Reads

### Strategy
Assemble unmapped reads de novo, then BLAST/diamond against viral protein database. No reference needed.

### Steps
1. Extract unmapped reads (same as D1)
2. Assemble: `spades.py --rnaviral -1 R1.fq.gz -2 R2.fq.gz -o spades_viral/`  (or Trinity)
3. BLAST contigs: `diamond blastx -d viral_proteins.dmnd -q contigs.fasta -o viral_hits.tsv`

### Recommendation
**Not recommended for OSCAR** — high compute cost, low yield for short-read scRNA-seq where viral coverage is sparse. Reserve for targeted investigation of specific samples with strong viral signal.

---

## Recommended Implementation Order

1. **Method A** — lowest effort, biggest payoff for known human-tropic herpesviruses common in immune cell biology (EBV, CMV, HHV-6/7/8). Build once, use forever.
2. **Method C** — add Kraken2 as a screening step in the MultiQC report. Runs fast, integrates with existing MultiQC call, requires no changes to counting logic.
3. **Method B or D** — only if per-cell viral resolution is needed beyond Method A, or for novel/unexpected viruses.

---

## Phase 1: Viral Reference Database Construction (Methods A + C)

### Tasks
1. Write `bash/build_viral_reference.sh`: download selected viral genomes, clean FASTA headers, soft-mask, concatenate to `hg38_viral.fa` + `hg38_viral.gtf`, run `cellranger mkref`
2. Write `bash/build_kraken2_db.sh`: build standard Kraken2 viral database
3. Document database versions and download dates in `assets/viral_references/README.md`
4. Add `params.ref_human_viral` and `params.kraken2_viral_db` to `nextflow.config` (default: null = disabled)

### Verification
- `cellranger mkref` completes without error
- `grep "^virus_" hg38_viral.gtf | wc -l` > 0
- `kraken2 --db kraken2_viral_db --version` exits 0

---

## Phase 2: Pipeline Integration

### Tasks
1. **Method A** — in `subworkflows/count_gex.nf`: replace `params.ref_human` with `params.ref_human_viral ?: params.ref_human` in the multi_config GEX reference line
2. **Method C** — new process `KRAKEN2_VIRAL` in `modules/qc.nf`; add to GEX ch_fastqs path in `main.nf`; add reports to MULTIQC inputs

### New params
```groovy
// ── Viral detection ────────────────────────────────────────────────────────
ref_human_viral    = null   // path to hg38 + viral combined cellranger reference
kraken2_viral_db   = null   // path to Kraken2 viral database
```

### Samplesheet: no changes required

---

## Phase 3: Verification

1. Test Method A: run `--from_cellranger` on an EBV-transformed B cell line sample; confirm `virus_EBV_*` genes appear in feature matrix
2. Test Method C: run Kraken2 on a known CMV+ sample; confirm CMV reads classified; MultiQC report shows viral section
3. Check no regression in non-viral samples (no spurious viral reads in healthy T cell controls)

---

## Anti-patterns to Avoid

- Do not add viral genomes to ATAC references — viral detection in ATAC is not meaningful (chromatin accessibility, not transcription)
- Do not alias viral contig names without `virus_` prefix — collision with human chromosome names breaks cellranger
- Do not run de novo assembly (Method E) at scale — prohibitively expensive for standard QC
- Do not use Kraken2 confidence < 0.05 — too many false positives from human reads with k-mer similarity to viral sequences
- Do not include endogenous retroelements (HERV) in the viral detection set without explicit intent — they are human genome features, not exogenous viral infections
