# OSCAR

**Ollie's Single Cell Analysis for the Romagnani Lab**

OSCAR takes a raw BCL folder and a samplesheet, and gives you demultiplexed
FASTQs, count matrices, and QC output. It runs on SLURM through Nextflow DSL2,
with every tool in an Apptainer container.

Full documentation: [ollieeknight.github.io/OSCAR](https://ollieeknight.github.io/OSCAR/)

## Supported assays

| Assay | Modalities | Counted with |
|-------|-----------|--------------|
| GEX | GEX, VDJ-T, VDJ-B, CRISPR | cellranger multi |
| CITE | GEX, ADT, HTO | cellranger multi |
| Flex | GEX | cellranger multi, cyto, or both |
| ATAC | ATAC | cellranger-atac |
| Multiome | GEX, ATAC | cellranger multi, cellranger-atac |
| DOGMA | GEX, ATAC, ADT, HTO | cellranger multi, cellranger-atac |
| ASAP | ATAC, ADT, HTO | cellranger-atac, kallisto/bustools |

Beyond counting, OSCAR runs falco and MultiQC on reads, cellbender and scrublet
on GEX, AMULET and mgatk2 and MACS3 on ATAC, and cellsnp-lite with vireo when a
library has more than one donor. Viral detection and RNA velocity are off by
default; enable them with `--extras viral,velocity`.

## Getting started

```bash
git clone https://github.com/ollieeknight/OSCAR
cd OSCAR
```

Build references with the scripts in `assets/references/`, then point
`nextflow.config` at your cluster paths and SLURM partitions. See
[Installation](https://ollieeknight.github.io/OSCAR/guide/installation/).

Run one flowcell:

```bash
nextflow run main.nf -profile slurm \
    --samplesheet   /path/to/metadata.csv \
    --bcl_dir       /path/to/R463_bcl \
    --adt_files_dir /path/to/adt_files \
    --outdir        /path/to/results \
    --run_name      R463
```

Merge libraries sequenced across several flowcells with `--extra_bcl_dirs` and
`--extra_samplesheets`; each flowcell's FASTQs are written beside its own BCL
folder. Skip demultiplexing with `--run_from fastq`, or run QC alone with
`--run_from cellranger`.

Build a samplesheet with the
[metadata generator](https://ollieeknight.github.io/OSCAR/tools/metadata_generator.html),
or copy `assets/example_metadata.csv`.

## Checking your setup

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
nextflow run tests/test_lib.nf
```

The first parses every file. The second exercises index kit loading, chemistry
lookup, samplesheet parsing, and cellranger config generation against the
shipped example. Neither needs a cluster.

## Repository layout

| Path | Contents |
|------|----------|
| `main.nf` | Routing and channel logic |
| `nextflow.config` | Parameters, profiles, per-process resources |
| `lib/` | Pure functions: chemistry registry, index kits, samplesheet parsing |
| `modules/` | One process per tool |
| `subworkflows/` | Process chains |
| `tests/` | Self-checks that run without a cluster |
| `assets/` | Index kit CSVs, reference build scripts, example samplesheets |
| `docs-src/` | Documentation source |
| `docs/` | Built site (generated, see below) |
| `original/` | Legacy bash implementation, reference only |

## Building the docs

GitHub Actions builds and deploys the site on every push to `main` that touches
`docs-src/`. Pull requests build without deploying, so a broken link fails the
check rather than the live site.

Preview locally:

```bash
pip install -r docs-src/requirements.txt
mkdocs serve                 # localhost:8000
```

## Contact

`oliver.knight@charite.de`
