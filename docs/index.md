# Ollie's Single Cell Analysis for the Romagnani lab (OSCAR)

<div style="text-align: center;">
    <img src="assets/images/oscar.jpg" alt="OSCAR - The best boy" width="300" style="border-radius: 8px; margin: 20px 0;">
</div>

OSCAR is a comprehensive pipeline for single-cell sequencing analysis, developed for the Romagnani lab at Charité. It provides tools and workflows for processing CITE-seq, DOGMA-seq, ASAP-seq, and Tapestri experiments.

## Quick Start

### 1. Clone the Repository

```bash
cd $HOME/work/bin/
git clone https://github.com/ollieeknight/OSCAR
```

### 2. Download Apptainer Images

There are two Apptainer images:

- `oscar-count.sif` - For counting steps
- `oscar-qc.sif` - For post-counting quality control

```bash
mkdir -p ${TMPDIR}/OSCAR
apptainer pull library://romagnanilab/oscar/oscar-count:latest --dir ${TMPDIR}/OSCAR/
apptainer pull library://romagnanilab/oscar/oscar-qc:latest --dir ${TMPDIR}/OSCAR/
```

### 3. Reference Genomes

Reference genomes are also required. Build steps can be found under the `reference/` directory.

## Features

- **Automated Metadata Generation** - Create standardized metadata files for your experiments
- **Feature Barcode References** - Generate barcode references for Cell Ranger, Kallisto, and Tapestri
- **Streamlined Processing** - Bash scripts to automate the entire workflow
- **Quality Control** - Comprehensive QC steps built into the pipeline

## Supported Assays

- **CITE-seq** - Cellular Indexing of Transcriptomes and Epitopes by Sequencing
- **DOGMA-seq** - DNA, Oligonucleotide, Gene expression, and Methylation Analysis
- **ASAP-seq** - ATAC with Select Antigen Profiling by sequencing
- **Tapestri** - Single-cell DNA and protein sequencing

## Getting Help

For questions, suggestions, or issues:

- **Email**: [oliver.knight@charite.de](mailto:oliver.knight@charite.de)
- **GitHub Issues**: [Report a bug or request a feature](https://github.com/ollieeknight/OSCAR/issues)

## Next Steps

- [Getting Started](getting-started.md) - Detailed installation and setup guide
- [Metadata Generator](tools/metadata-generator.html) - Create metadata files for your experiments
- [Feature Barcode Generator](tools/adt-generator.html) - Generate feature barcode references
- [Functions Reference](reference/functions.md) - Detailed explanation of all OSCAR functions
