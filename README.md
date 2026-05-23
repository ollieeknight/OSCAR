# OSCAR

**Ollie's Single Cell Analysis for the Romagnani Lab**

OSCAR is a comprehensive pipeline designed for processing single-cell RNA, ATAC, and multiome sequencing data.

## Features
- **Comprehensive metadata tracking**: Generate your metadata file [here](https://ollieeknight.github.io/OSCAR/) and keep track of sequencing runs.
- **FASTQ demultiplexing and QC**: Takes a raw bcl folder and demultiplexes FASTQ files, performs `falco` and `multiqc`.
- **scRNA pipeline**: Utilises `cellranger` for counting (with and without ADT/HTO), performs ambient-RNA correction with `CellBender`, and donor genotyping demultiplexing with `cellsnp-lite` and `vireo`.
- **scATAC pipeiline**: Utilises `cellranger` for counting, AMULET for doublet prediction, donor genotyping demultiplexing with `cellsnp-lite` and `vireo`, and for ASAP-seq performs ADT/HTO counting with `kallisto`.
- **Multiome pipeline**: All of the above
- **Multiple sequencing run integration**: Integrate sequenced libraries from several runs into one output folder.

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/ollieeknight/OSCAR
   cd OSCAR
   ```

2. Build reference genomes (see `reference/` for human and mouse build scripts).

3. Run the pipeline (BCL → FASTQ → count → QC):
   ```bash
   nextflow run main.nf -profile slurm \
       --samplesheet /path/to/metadata.csv \
       --bcl_dir     /path/to/BCL_folder \
       --outdir      results \
       --adt_files_dir /path/to/adt_csvs/
   ```

See `CLAUDE.md` for full parameter reference, samplesheet format, and invocation examples.

For any questions, please e-mail `oliver.knight@charite.de`.

> **Note:** Legacy bash scripts are archived in `old/bash/` for reference only.
