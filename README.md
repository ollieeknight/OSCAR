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

Beyond counting, OSCAR runs fastp and MultiQC on reads, cellbender and scrublet
on GEX, AMULET and mgatk2 and MACS3 on ATAC, and cellsnp-lite with vireo when a
library has more than one donor. Viral detection and RNA velocity are off by
default; enable them with `--extras viral,velocity`.

## Requirements

- Nextflow 24.04 or newer, and Apptainer
- A SLURM cluster (or run locally with `-profile standard`)
- Reference genomes, built with the scripts in `assets/references/`

Every tool runs in a container, so nothing else needs installing.

## Setting up

```bash
git clone https://github.com/ollieeknight/OSCAR
cd OSCAR
```

**Before your first run, edit `nextflow.config`.** The reference, whitelist and
scratch paths shipped in it point at the Romagnani Lab's filesystem, and none of
them will resolve anywhere else. Set `ref_human`, `ref_mouse`, the VDJ and SNP
references, `workDir`, and the SLURM partitions to your own. See
[Installation](https://ollieeknight.github.io/OSCAR/guide/installation/).

Build a samplesheet with the
[metadata generator](https://ollieeknight.github.io/OSCAR/tools/metadata_generator.html),
or copy `assets/example_metadata.csv`.

## Running

Every run needs `--samplesheet` and one input: `--bcl_dir` to demultiplex from
BCL (the default), `--fastq_dir` to start from FASTQs, or `--outs_dir` to run QC
on existing cellranger output.

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir     /path/to/R463_bcl
```

Output lands beside the input flowcell: FASTQs and read-level QC in
`{run}_fastq`, results in `{run}_outs`, both siblings of `{run}_bcl`.

### Common options

| Option | When you need it |
|--------|------------------|
| `--adt_files_dir` | CITE, ASAP or DOGMA runs, pointing at your feature barcode CSVs |
| `--run_name` | Naming the output something other than the BCL folder |
| `--outdir` | Writing results somewhere other than beside the flowcell |
| `--extra_bcl_dirs`, `--extra_samplesheets` | One library sequenced across several flowcells |
| `--run_from fastq\|cellranger` | Skipping demultiplexing, or running QC alone |
| `--extras viral,velocity` | Viral detection and RNA velocity, both off by default |
| `--flex_backend` | Flex runs: `cellranger`, `cyto` or `both` |

Merging flowcells writes each one's FASTQs beside its own BCL folder, and
collects the counts under a single `{run}_outs` named by `--run_name`.

Every parameter is listed in the
[parameter reference](https://ollieeknight.github.io/OSCAR/reference/parameters/).

## Checking your setup

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
bash tests/run_all.sh
```

The first parses every file. The second runs every self-check: index kit
loading, chemistry lookup, samplesheet parsing, OverrideCycles generation,
FASTQ staging and QC routing, demux QC, and the publish layout. Neither needs a
cluster, containers, or sequencing data.

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

## Citation and licence

OSCAR is released under the MIT licence; see [LICENSE](LICENSE). If you use it
in published work, cite it using the metadata in [CITATION.cff](CITATION.cff).

## Contact

Oliver Knight, `oliver.knight@charite.de`
