# OSCAR Agent Guidelines

Nextflow DSL2 single-cell pipeline. For Romagnani Lab (BIH Charité). Process scRNA-seq, scATAC-seq, Multiome, DOGMA, ASAP-seq, CITE-seq, Flex (Fixed RNA Profiling), VDJ, CRISPR. Run on SLURM, Apptainer container.

---

## Directory Structure

```
.
├── main.nf                    ← Pipeline entrypoint, channel routing, CLI params, preflight
├── nextflow.config            ← Parameters, containers, SLURM profile, resource allocations
├── lib/
│   ├── chemistry.nf           ← Chemistry registry (single source of truth for chemistries/whitelists)
│   ├── indexes.nf             ← 10x SI/DI kit index parsing & sequencer auto-detection
│   ├── samplesheet.nf         ← Samplesheet parsing, validation, ADT CSV resolution
│   └── multi_config.nf        ← cellranger multi config header & flex sample builders
├── modules/
│   ├── demux.nf               ← GENERATE_SAMPLESHEET, BCLCONVERT, VALIDATE_FASTQ, FALCO, MULTIQC
│   ├── count_gex.nf           ← CELLRANGER_MULTI (staged FASTQs, Python read filter, config gen)
│   ├── count_gex_cyto.nf      ← CYTO_FLEX, CYTO_RENAME_SAMPLES (Flex probe-level QC)
│   ├── flex_probe_convert.nf  ← FLEX_PROBE_PREPARE, FLEX_SAMPLE_PREPARE, barcode extraction
│   ├── count_atac.nf          ← CELLRANGER_ATAC (staged FASTQs, Python read filter)
│   ├── count_adt.nf           ← ASAP ADT kallisto/bustools pipeline
│   ├── qc_gex.nf              ← CELLBENDER (GPU), SCRUBLET
│   ├── qc_atac.nf             ← AMULET, MGATK2, MACS3
│   ├── genotype.nf            ← CELLSNP_LITE, VIREO (shared GEX/ATAC donor demux)
│   └── quant_extra.nf         ← VIRAL_DETECT, SIMPLEAF_VELOCITY
├── subworkflows/
│   ├── demux.nf               ← DEMUX: samplesheet gen, per-lane BCL Convert, FASTQ QC
│   ├── fastq_qc.nf            ← FASTQ_QC: pigz validation + Falco on R-reads
│   ├── count_gex.nf           ← COUNT_GEX: probe preparation, cellranger multi, optional cyto
│   ├── count_atac.nf          ← COUNT_ATAC: CELLRANGER_ATAC execution
│   ├── count_adt.nf           ← COUNT_ADT: ASAP kallisto/bustools processing
│   ├── qc_gex.nf              ← QC_GEX: CellBender → Scrublet → Genotyping
│   ├── qc_atac.nf             ← QC_ATAC: AMULET + mgatk2 + MACS3 → Genotyping
│   └── genotype.nf            ← GENOTYPE: cellsnp-lite → vireo wrapper
└── assets/
    ├── indexes/               ← 10x SI/DI kit CSV index definitions
    └── multiqc_config.yaml    ← MultiQC report configuration
```

---

## Pipeline Execution & Entry Points

### 1. BCL Demux Mode (Default)
```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir /path/to/BCL_dir \
    --outdir results \
    --run_name R463
```
- Many flowcell: use `--extra_bcl_dirs /path/to/R464` (plus `--extra_samplesheets` if sheets differ).
- FASTQ merge by `meta.library_id`, downstream, before count.

### 2. FASTQ Mode
```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from_fastq true \
    --fastq_dir /path/to/fastqs \
    --outdir results
```

### 3. CellRanger Outs Mode (QC Only)
```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from_cellranger true \
    --outs_dir /path/to/results \
    --outdir results
```

### Execution Limits
- `--run_until FASTQ`: demux + fastq QC only.
- `--run_until cellranger`: count only, skip QC.

---

## Core Routing Rules

| Assay | Modality | Route | Notes |
|---|---|---|---|
| GEX, CITE, DOGMA, Multiome, Flex | GEX, ADT, HTO, VDJ-T, VDJ-B, CRISPR | COUNT_GEX → QC_GEX | DOGMA ADT run cellranger multi, never kallisto |
| ATAC, DOGMA, Multiome, ASAP | ATAC | COUNT_ATAC → QC_ATAC | Standalone or paired ATAC count |
| ASAP | ADT, HTO | COUNT_ADT (kallisto) | Trigger only after ATAC cellranger done |
| Any | GENO | Skipped | Info only |

- **DOGMA Experiments**: GEX and ATAC run in separate Nextflow runs, integrate downstream in R.
- **Flex Experiments**: Use `flex_backend` (`cellranger`, `cyto`, or `both`). Need `--flex_probe_set` and `--flex_samples_file` when multiplexed.
- **Donor Genotyping**: `CELLSNP_LITE` → `VIREO` run only when `meta.n_donors > 1 && meta.species == 'human'`. Mouse sample always set `n_donors = 1`.

---

## Key Technical Mechanics

- **FASTQ Staging & Filtering**: FASTQ file stage via `path(fastqs, stageAs: "fastqs/{mod}/run_???/*")`. Python inside `CELLRANGER_MULTI` / `CELLRANGER_ATAC` filter empty placeholder run (`<10000` reads on R1), rename lane in order to `L001`, `L002`, etc.
- **Nextflow ArrayBag Materialization**: `groupTuple` emit `nextflow.util.ArrayBag`. No do collection op inside process script block. Materialize inside subworkflow `.map {}` block:
  ```groovy
  .map { key, metas ->
      def ml = []
      metas.each { m -> ml << m }
      [key, ml]
  }
  ```
- **Sequencer Auto-Detection**: `detect_sequencer` look at `<Instrument>` in `RunInfo.xml`, per BCL folder:
  - `VH`, `NDX`, `FS`, `LH` → `novaseq_x` (i5 forward).
  - `A`, `NB`, `NS`, `MN` → `novaseq6000` (i5 reverse complement).
- **Chemistry Registry**: `lib/chemistry.nf` = one truth source for chemistry, family, barcode whitelist, simpleaf chemistry.
- **ADT CSV Resolution**: Look for `{samplesheet_dir}/adt_files/{adt_file}.csv`, then `{samplesheet_dir}/../adt_files/`, then `--adt_files_dir`.
- **Optional Analyses**:
  - `VIRAL_DETECT`: simpleaf viral quant, on unassigned BAM read (gate on `params.viral_piscem_index`).
  - `SIMPLEAF_VELOCITY`: simpleaf spliced/unspliced quant on GEX, use CellBender filtered barcode (gate on `params.run_velocity`, skip for Flex).

---

## Skills to Use

In this repo, follow these agent skill rule:

- **`investigate-first`**: Use before touch code, when diagnose failed pipeline run, Slurm job crash, missing input file, or bad channel tuple. Read task `.command.log`, `.command.err`, and `.command.sh`.
- **`surgical-patch`**: Use for bug fix in module, library script, or config file. Keep edit small, focus on broke part, no touch surround logic.
- **`safe-refactor`**: Use when change subworkflow topology, channel operator, or process signature. Check emit tuple match downstream input need.
- **`ponytail` / `ponytail-review`**: Enforce minimal, YAGNI. Favor native Nextflow DSL2 operator and Python stdlib over extra dependency or wrapper abstraction.
- **`cavecrew`**: Send deep investigation, multi-module search, or big log audit to investigator subagent. Keep context clean.
- **`verify-and-stop`**: Test syntax (`nextflow config` / dry-run), check acceptance criteria, no speculative change.