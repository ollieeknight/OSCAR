# Quickstart

## 1. Prepare files

Before starting OSCAR, prepare:

1. A completed BCL folder.
2. One metadata CSV with the [metadata generator](../tools/metadata_generator.html).
3. For ADT or HTO libraries, a reference CSV with the [feature barcode generator](../tools/adt_generator.html). Save it in `adt_files/` beside `metadata.csv`.

## 2. Run OSCAR

```bash
nextflow run main.nf -profile slurm \
  --samplesheet /path/to/metadata.csv \
  --bcl_dir /path/to/R463_bcl \
  --outdir /path/to/results \
  --run_name R463
```

Add `--adt_files_dir /path/to/adt_files` only when references are stored elsewhere. Omit `--run_name` to use the BCL folder name.

## 3. Find results

```
{outdir}/
├── {run_name}_fastq/multiqc/multiqc_report.html
└── {run_name}_outs/
    ├── {library_id}/outs/        GEX, CITE, Flex, Multiome, DOGMA
    └── {library_id}_ATAC/outs/   ATAC, Multiome, DOGMA, ASAP
```

Open MultiQC first. Each library folder has counts and relevant QC output.

## Common variations

**Same libraries across several flowcells**

```bash
--extra_bcl_dirs /path/to/R462_bcl,/path/to/R463_bcl \
--extra_samplesheets /path/to/R462.csv,/path/to/R463.csv
```

Keep both lists in the same order. Omit `--extra_samplesheets` when every flowcell uses the same metadata CSV.

**Start from FASTQs**

```bash
--from_fastq --fastq_dir /path/to/fastqs
```

**Stop after a stage**

```bash
--run_until FASTQ
--run_until cellranger
```

Use underscores in parameter names. See [Parameters](../reference/parameters.md) for uncommon options.
