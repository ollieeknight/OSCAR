# OSCAR Agent Guidelines

Nextflow DSL2 single-cell pipeline for the Romagnani Lab (BIH Charité). Processes scRNA-seq, scATAC-seq, Multiome, DOGMA, ASAP-seq, CITE-seq, Flex (Fixed RNA Profiling), VDJ and CRISPR data. Runs on SLURM under Apptainer.

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
│   ├── demux.nf               ← GENERATE_SAMPLESHEET, CLEAN_FASTQ_DIR, BCLCONVERT, DEMUX_QC,
│   │                            CELLRANGER_MQC, FASTP, MULTIQC, VALIDATE_FASTQ
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
│   ├── fastq_qc.nf            ← FASTQ_QC: fastp on R-reads, pigz on everything else
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
- Several flowcells: add `--extra_bcl_dirs /path/to/R464`, plus `--extra_samplesheets` when the sheets differ.
- FASTQs publish beside their source flowcell at `{bcl_dir.parent}/{run}_fastq`, not under `--outdir`.
- FASTQs merge by `meta.library_id` downstream, before counting.

### 2. FASTQ Mode
```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --run_from fastq \
    --fastq_dir /path/to/fastqs \
    --outdir results
```

### 3. CellRanger Outs Mode (QC Only)
```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --run_from cellranger \
    --outs_dir /path/to/results \
    --outdir results
```

### Execution Limits
- `--run_from bcl` (default), `fastq` or `cellranger` selects the entry point.
- `--extras velocity,viral` opts in to supplementary quantification.

---

## Core Routing Rules

| Assay | Modality | Route | Notes |
|---|---|---|---|
| GEX, CITE, DOGMA, Multiome, Flex | GEX, ADT, HTO, VDJ-T, VDJ-B, CRISPR | COUNT_GEX → QC_GEX | DOGMA ADT runs through cellranger multi, never kallisto |
| ATAC, DOGMA, Multiome, ASAP | ATAC | COUNT_ATAC → QC_ATAC | Standalone or paired ATAC counting |
| ASAP | ADT, HTO | COUNT_ADT (kallisto) | Triggered once ATAC cellranger finishes |
| Any | GENO | Skipped | Informational only |

- **DOGMA experiments**: GEX and ATAC run as separate Nextflow runs and are integrated downstream in R.
- **Flex experiments**: set `flex_backend` to `cellranger`, `cyto` or `both`. Multiplexed runs need `--flex_probe_set` and `--flex_samples_file`.
- **Donor genotyping**: `CELLSNP_LITE` and `VIREO` run only when `meta.n_donors > 1 && meta.species == 'human'`. Mouse samples always set `n_donors = 1`.

---

## Key Technical Mechanics

- **FASTQ staging and filtering**: files stage via `path(fastqs, stageAs: "fastqs/{mod}/run_???/*")`, which gives each file its own `run_???` dir. `bin/stage_fastqs.py` pairs reads on the FASTQ header, since two flowcells can stage identically-named files and the dir name records list position rather than flowcell. `CELLRANGER_MULTI` and `CELLRANGER_ATAC` then drop placeholder runs (under 10000 reads on R1) and renumber the survivors `L001`, `L002` and so on.
- **ArrayBag materialisation**: `groupTuple` emits a `nextflow.util.ArrayBag`. Do not run collection operations on one inside a process script block; convert it in the subworkflow `.map {}` first:
  ```groovy
  .map { key, metas -> [key, metas.collect()] }
  ```
- **Sequencer auto-detection**: `detect_sequencer` reads `<Instrument>` from `RunInfo.xml`, once per BCL folder:
  - `VH`, `NDX`, `FS`, `LH` → `novaseq_x` (i5 forward).
  - `A`, `NB`, `NS`, `MN` → `novaseq6000` (i5 reverse complement).
- **Chemistry registry**: `lib/chemistry.nf` is the single source of truth for chemistry, family, barcode whitelist and simpleaf chemistry.
- **ADT CSV resolution**: checks `{samplesheet_dir}/adt_files/{adt_file}.csv`, then `{samplesheet_dir}/../adt_files/`, then `--adt_files_dir`.
- **Optional analyses**:
  - `VIRAL_DETECT`: simpleaf viral quantification over unassigned BAM reads, gated on `--extras viral`. `bamtofastq_bin` may be an https URL, which Nextflow stages and caches.
  - `SIMPLEAF_VELOCITY`: simpleaf spliced/unspliced quantification on GEX using the CellBender barcode list, gated on `--extras velocity` and skipped for Flex.

---

## Working in This Repo

- Diagnose a failed run from the task directory first: read `.command.log`,
  `.command.err` and `.command.sh` before changing code.
- Keep bug fixes small and local to the broken part.
- Changing subworkflow topology, a channel operator or a process signature
  means re-checking that the emitted tuple still matches what downstream
  consumes.
- Prefer native DSL2 operators and the Python standard library over new
  dependencies or wrapper layers.
- Verify before claiming done: `bash tests/run_all.sh` runs every self-check
  without a cluster, containers or sequencing data.
