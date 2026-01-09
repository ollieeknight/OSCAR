# Getting Started

This guide will walk you through setting up OSCAR for single-cell sequencing analysis.

## Prerequisites

- Access to a Linux/Unix system (BIH cluster recommended)
- Apptainer/Singularity installed
- Basic command-line knowledge

## Installation

### 1. Clone the Repository

```bash
cd $HOME/work/bin/
git clone https://github.com/ollieeknight/OSCAR
cd OSCAR
```

### 2. Download Container Images

OSCAR uses two Apptainer containers for different stages of the analysis:

#### oscar-count.sif
Used for alignment and counting:
```bash
mkdir -p ${TMPDIR}/OSCAR
apptainer pull library://romagnanilab/oscar/oscar-count:latest --dir ${TMPDIR}/OSCAR/
```

#### oscar-qc.sif
Used for quality control and downstream analysis:
```bash
apptainer pull library://romagnanilab/oscar/oscar-qc:latest --dir ${TMPDIR}/OSCAR/
```

### 3. Set Up Reference Genomes

Reference genomes are required for alignment. Build scripts are provided in the `reference/` directory:

#### Human Reference (GRCh38)
```bash
cd reference/
bash human.sh
```

#### Mouse Reference (GRCm39)
```bash
cd reference/
bash mouse.sh
```

#### Kallisto Index (for ASAP-seq)
```bash
cd reference/
bash human_kallisto.sh
```

## Project Setup

### 1. Create Project Directory

```bash
PROJECT_ID="my_experiment"
DIR_PREFIX="${HOME}/scratch/ngs"
mkdir -p ${DIR_PREFIX}/${PROJECT_ID}
cd ${DIR_PREFIX}/${PROJECT_ID}
```

### 2. Prepare Metadata File

Create a `metadata.csv` file or use the [Metadata Generator](tools/metadata-generator.html) tool.

Example metadata structure:
```csv
assay,experiment_id,historical_number,replicate,modality,chemistry,index_type,index,species,n_donors,adt_file
CITE,EXP001,H001,R1,RNA,3prime,dual,SI-TT-A1,Human,1,adt_totalseq_a.csv
```

### 3. Add FASTQ Files

Place your FASTQ files in:
```bash
${DIR_PREFIX}/${PROJECT_ID}/fastq/
```

## Running the Pipeline

### Step 1: Process Metadata
```bash
bash ${HOME}/work/bin/OSCAR/bash/01_process_metadata.sh \
    --project-id ${PROJECT_ID} \
    --dir-prefix ${DIR_PREFIX}
```

### Step 2: Process FASTQ Files
```bash
bash ${HOME}/work/bin/OSCAR/bash/02_fastq.sh \
    --project-id ${PROJECT_ID} \
    --dir-prefix ${DIR_PREFIX}
```

### Step 3: Process Libraries
```bash
bash ${HOME}/work/bin/OSCAR/bash/03_process_libraries.sh \
    --project-id ${PROJECT_ID} \
    --dir-prefix ${DIR_PREFIX}
```

### Step 4: Count
```bash
bash ${HOME}/work/bin/OSCAR/bash/04_count.sh \
    --project-id ${PROJECT_ID} \
    --dir-prefix ${DIR_PREFIX}
```

### Step 5: Quality Control
```bash
bash ${HOME}/work/bin/OSCAR/bash/05_quality_control.sh \
    --project-id ${PROJECT_ID} \
    --dir-prefix ${DIR_PREFIX}
```

## Configuration

### Custom Settings

Edit `bash/config.sh` to customize:

- Reference genome paths
- Container image locations
- Default parameters

### Environment Variables

```bash
export OSCAR_HOME="${HOME}/work/bin/OSCAR"
export OSCAR_IMAGES="${TMPDIR}/OSCAR"
```

## Troubleshooting

### Common Issues

!!! warning "FASTQ files not found"
    Ensure FASTQ files follow the naming convention expected by your sequencing platform.

!!! warning "Memory errors"
    Increase memory allocation in your job submission script.

!!! tip "Checking logs"
    Log files are created in `${DIR_PREFIX}/${PROJECT_ID}/logs/`

## Next Steps

- [Metadata Generator](tools/metadata-generator.html) - Create properly formatted metadata
- [Feature Barcode Generator](tools/adt-generator.html) - Generate ADT/HTO reference files
- [Functions Reference](reference/functions.md) - Detailed script documentation

## Support

For help or questions:

- **Email**: [oliver.knight@charite.de](mailto:oliver.knight@charite.de)
- **GitHub Issues**: [Submit an issue](https://github.com/ollieeknight/OSCAR/issues)
