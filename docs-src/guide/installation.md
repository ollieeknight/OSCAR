# Installation

Do this once per cluster. Lab users normally only need a working OSCAR directory and `slurm` profile.

## Requirements

OSCAR needs Nextflow, SLURM, and Apptainer on compute nodes. OSCAR downloads other software automatically.

```bash
nextflow -version
apptainer --version
```

## Get OSCAR

```bash
git clone https://github.com/ollieeknight/OSCAR
cd OSCAR
```

## Configure your cluster

Set these values in `nextflow.config` before the first run:

- `workDir`: fast scratch space.
- `apptainer.cacheDir`: image cache.
- `apptainer.runOptions`: folders OSCAR can read and write.
- `process.queue`: SLURM queues, including the GPU queue for CellBender.
- `ref_*`, `snp_vcf`, `atac_whitelist`: reference locations.

Build the two local containers once. Set `container_cyto` and `container_cellbender` to the resulting files.

```bash
apptainer build /path/to/cache/cyto.sif containers/cyto.def
```

## Check setup

```bash
nextflow run tests/test_lib.nf
```

Expected result: `OK: all lib/ self-checks passed`.
