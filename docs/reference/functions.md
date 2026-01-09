# OSCAR Functions Reference

This page provides detailed documentation for all OSCAR bash scripts.

## 01_process_metadata.sh

**Description**: This script processes the metadata file for a given project. It reads the metadata file, determines the run type, and creates index files for each library based on the metadata.

**Inputs**:

- `--project-id`: The project ID
- `--dir-prefix`: The directory prefix (default: `${HOME}/scratch/ngs`)
- `--metadata-file-name`: The metadata file name (default: `metadata.csv`)

**Outputs**: Index files for each library in the `indices` folder

**Example**:
```bash
bash 01_process_metadata.sh \
    --project-id MY_PROJECT \
    --dir-prefix ${HOME}/scratch/ngs \
    --metadata-file-name metadata.csv
```

---

## 02_fastq.sh

**Description**: This script performs FASTQ demultiplexing for a given project. It uses the index files created by `01_process_metadata.sh` to run the cellranger command and generate FASTQ files.

**Inputs**:

- `--project-id`: The project ID
- `--dir-prefix`: The directory prefix (default: `${HOME}/scratch/ngs`)
- `--metadata-file-name`: The metadata file name (default: `metadata.csv`)

**Outputs**: FASTQ files and quality control reports in the `fastq` folder

**Example**:
```bash
bash 02_fastq.sh \
    --project-id MY_PROJECT \
    --dir-prefix ${HOME}/scratch/ngs
```

---

## 03_process_libraries.sh

**Description**: This script processes the libraries for a given project. It reads the metadata file, determines the run type, and handles different modes (GEX or ATAC) to prepare the libraries for counting.

**Inputs**:

- `--project-id`: The project ID
- `--dir-prefix`: The directory prefix (default: `${HOME}/scratch/ngs`)
- `--gene-expression-options`: Options for gene expression processing (comma-separated)
- `--vdj-options`: Options for VDJ processing (comma-separated)
- `--adt-options`: Options for ADT processing (comma-separated)
- `--metadata-file-name`: The metadata file name (default: `metadata.csv`)

**Outputs**: Prepared libraries for counting in the `libraries` folder

**Example**:
```bash
bash 03_process_libraries.sh \
    --project-id MY_PROJECT \
    --dir-prefix ${HOME}/scratch/ngs \
    --gene-expression-options "expect-cells=10000"
```

---

## 04_count.sh

**Description**: This script submits counting jobs for each library in a given project. It reads the library files, determines the run type, and submits the appropriate counting jobs to SLURM.

**Inputs**:

- `--project-id`: The project ID
- `--dir-prefix`: The directory prefix (default: `${HOME}/scratch/ngs`)
- `--metadata-file-name`: The metadata file name (default: `metadata.csv`)

**Outputs**: Counting job submissions and logs in the `outs` folder

**Example**:
```bash
bash 04_count.sh \
    --project-id MY_PROJECT \
    --dir-prefix ${HOME}/scratch/ngs
```

---

## 05_quality_control.sh

**Description**: This script performs quality control for a given project. It checks the necessary OSCAR containers, submits jobs for ambient RNA removal, genotyping, and other quality control tasks.

**Inputs**:

- `--project-id`: The project ID
- `--dir-prefix`: The directory prefix (default: `${HOME}/scratch/ngs`)
- `--metadata-file-name`: The metadata file name (default: `metadata.csv`)

**Outputs**: Quality control job submissions and logs in the `outs` folder

**Example**:
```bash
bash 05_quality_control.sh \
    --project-id MY_PROJECT \
    --dir-prefix ${HOME}/scratch/ngs
```

---

## Common Options

All scripts support the following common options:

| Option | Description | Default |
|--------|-------------|---------|
| `--project-id` | Unique identifier for your project | *Required* |
| `--dir-prefix` | Base directory for analysis | `${HOME}/scratch/ngs` |
| `--metadata-file-name` | Name of metadata CSV file | `metadata.csv` |

## Workflow Order

The scripts should be run in numerical order:

1. **01_process_metadata.sh** - Parse metadata and create index files
2. **02_fastq.sh** - Demultiplex and generate FASTQ files  
3. **03_process_libraries.sh** - Prepare libraries for counting
4. **04_count.sh** - Run counting (Cell Ranger/Kallisto)
5. **05_quality_control.sh** - Perform QC and filtering

## Tips

!!! tip "Parallel Processing"
    Most scripts can handle multiple libraries automatically and submit parallel jobs to SLURM.

!!! warning "Check Logs"
    Always check log files in the `logs/` directory for any errors or warnings.

!!! info "Resuming Failed Jobs"
    If a job fails, you can typically re-run just that step without starting over.
