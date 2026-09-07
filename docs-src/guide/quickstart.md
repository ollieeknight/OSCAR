# Quickstart

## One flowcell, start to finish

```bash
nextflow run main.nf -profile slurm \
    --samplesheet   /path/to/metadata.csv \
    --bcl_dir       /path/to/R463_bcl \
    --adt_files_dir /path/to/adt_files \
    --outdir        /path/to/results \
    --run_name      R463
```

Drop `--run_name` and OSCAR derives it from the BCL folder name, stripping a
trailing `_bcl`. Drop `--adt_files_dir` when no library has ADT or HTO.

## Where output lands

```
{outdir}/
├── {run_name}_fastq/
│   ├── *.fastq.gz
│   ├── falco/
│   └── multiqc/multiqc_report.html
└── {run_name}_outs/
    ├── {library_id}/outs/           cellranger multi
    │   ├── cellbender/
    │   ├── scrublet/
    │   └── vireo/
    └── {library_id}_ATAC/outs/      cellranger-atac
        ├── AMULET/
        ├── mgatk2/
        ├── peaks/
        └── ADT/                     ASAP only
```

## Several flowcells, one set of libraries

Sequence the same libraries across four flowcells and you want one count matrix
per library, not four. Pass the extra BCL folders and their samplesheets:

```bash
nextflow run main.nf -profile slurm \
    --samplesheet ${BASE}/R502_scripts/metadata/metadata.csv \
    --bcl_dir     ${BASE}/R502_bcl \
    --run_name    R502 \
    --extra_bcl_dirs     /path/R462_bcl,/path/R463_bcl,/path/R501_bcl \
    --extra_samplesheets /path/R462.csv,/path/R463.csv,/path/R501.csv \
    --adt_files_dir      ${BASE}/R502_scripts/adt_files \
    --outdir             ${BASE}
```

!!! warning "Order matters"

    `--extra_bcl_dirs` and `--extra_samplesheets` must be the same length and
    the same order. OSCAR pairs them positionally and errors on a length
    mismatch.

OSCAR detects the sequencer for each BCL folder on its own, so you can mix a
NovaSeq X run with a NovaSeq 6000 run and the i5 orientation comes out right for
both.

Give the same samplesheet for every flowcell by omitting
`--extra_samplesheets`.

## Stop early

```bash
--run_until FASTQ       # demultiplex and QC, then stop
--run_until cellranger  # count, then stop before QC
```

## Start from existing files

=== "From FASTQs"

    ```bash
    nextflow run main.nf -profile slurm \
        --samplesheet /path/to/metadata.csv \
        --from_fastq --fastq_dir /path/to/fastqs \
        --outdir /path/to/results
    ```

    OSCAR matches FASTQs by filename prefix against the `id` it builds from
    each samplesheet row, so `CITE_PBMC_exp1_libA_GEX_S1_L001_R1_001.fastq.gz`
    matches the row that produces `CITE_PBMC_exp1_libA_GEX`.

=== "From cellranger output"

    ```bash
    nextflow run main.nf -profile slurm \
        --samplesheet /path/to/metadata.csv \
        --from_cellranger --outs_dir /path/to/{run_name}_outs \
        --outdir /path/to/results
    ```

    Runs QC alone against finished counts.

!!! warning "These conflict"

    `--from_fastq` and `--from_cellranger` cannot be combined, and
    `--run_until` does nothing alongside `--from_cellranger`.

## Optional steps

```bash
--run_velocity true                          # spliced/unspliced with simpleaf
--viral_piscem_index /path/to/piscem/index   # viral transcript detection
```

Velocity skips Flex libraries, since probe-based chemistry carries no intronic
signal. Viral detection runs on human libraries only.
