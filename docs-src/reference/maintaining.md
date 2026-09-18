# Maintaining OSCAR

For whoever runs and changes OSCAR. Lab members do not need this page.

## Where everything lives

OSCAR runs on the **Charité SC cluster**, not BIH-CUBI. Three storage areas
matter, each reached through a symlink in the home directory:

| Symlink | Real path | Holds |
|---|---|---|
| `~/work` | `/sc-projects/sc-proj-cc12-ag-romagnani` | Code, references, containers |
| `~/scratch` | `/sc-scratch/sc-scratch-cc12-ag-romagnani` | Work directories, image cache |
| `~/f-cc12-ag-romagnani` | `/charite-store-f/f-cc12-ag-romagnani` | Sequencing data |

The same home directory also has `~/genome` and `~/wes` pointing into the
sequencing store. Those belong to other pipelines.

`nextflow.config` writes these as real paths, not symlinks, so they resolve for
SLURM jobs regardless of whose account runs them. Taking OSCAR over means
getting group access to the storage, then changing `/home/knighto` in
`apptainer.runOptions` to your own home.

OSCAR itself lives at `~/work/bin/OSCAR`, alongside the other lab pipelines and
the tools they need: `nextflow`, `apptainer`, `cellranger-10.0.0` and
`cellranger-atac-2.2.0`. Run it from there rather than cloning your own copy,
so everyone runs the same commit. Counting itself happens in containers; the
local cellranger installs only supply `atac_whitelist` and
`tenx_barcodes_dir`.

Nextflow's work directory is `~/scratch/nf_work_oscar`, and Apptainer caches
images in `~/scratch/apptainer_cache`. Both are scratch and are not backed up.

References sit under `~/work/ref/hs` and `~/work/ref/mm`: genome indexes, the
1000 Genomes SNP VCF for donor demultiplexing, and the 10x Flex probe set. The
build scripts that produced them are in `assets/references/`.

## Configuring for a new account or cluster

Everything cluster-specific is in `nextflow.config`:

| Setting | What to change |
|---|---|
| `workDir` | Fast scratch, not backed up |
| `apptainer.cacheDir` | Where images are cached |
| `apptainer.runOptions` | Bind mounts. Every path OSCAR reads or writes needs one |
| `process.queue` | The normal SLURM partition |
| `ref_*`, `snp_vcf`, `atac_whitelist` | Reference locations |
| `container_*` | Images. Most pull automatically |

The `CELLBENDER` process overrides `queue` to `gpu` and needs `--gres=gpu:1`;
`VIRAL_DETECT` and `SIMPLEAF_VELOCITY` require AMD Genoa or Turin nodes.

Two containers are built locally rather than pulled, because no public image
exists. Build once and point `container_cyto` and `container_cellbender` at the
results:

```bash
apptainer build ~/scratch/apptainer_cache/cyto.sif containers/cyto.def
apptainer build ~/scratch/apptainer_cache/cellbender_sc1.sif containers/cellbender.def
```

The filenames must match `container_cyto` and `container_cellbender`.

## How a run flows

1. Parse the metadata CSV, validate it, expand index kit codes to sequences.
2. Demultiplex BCL to FASTQ with bcl-convert, one job per demux group per lane.
3. QC every FASTQ: fastp on reads, `pigz -t` on the rest, then MultiQC.
4. Count each library with the tool its assay needs.
5. Run whichever QC applies to the library's modalities.

Counting routes by assay: `cellranger multi` for GEX, CITE, Flex, Multiome and
DOGMA; `cellranger-atac` for ATAC and the ATAC half of Multiome, DOGMA and
ASAP; kallisto/bustools for ASAP antibody libraries.

## Repository map

| Path | Contents |
|---|---|
| `main.nf` | Routing and channel logic |
| `nextflow.config` | Parameters, profiles, per-process resources |
| `lib/` | Pure functions: chemistry, index kits, samplesheet parsing |
| `modules/` | One process per tool |
| `subworkflows/` | Process chains |
| `tests/` | Self-checks that run without a cluster |
| `assets/` | Index kit CSVs, reference build scripts, examples |
| `docs-src/` | This site's source |
| `original/` | Legacy bash implementation, reference only |

## Resources and retries

Every process has a `withName` block in the `slurm` profile. Processes that
can hit a ceiling scale memory and walltime by `task.attempt`, so a retry asks
for more than the attempt that died. Retries fire only on resource kills:
out-of-memory, walltime, and the signals a killed child raises. Anything else
fails immediately rather than burning a slot.

Single-threaded tools are given one CPU on purpose. Check before raising one:
`macs3`, `AMULET`, `bustools correct` and `bustools count` take no thread
argument, and `scrublet` has its BLAS threads capped to `task.cpus` in the
process itself.

## The generators

The two HTML tools are plain JavaScript with no build step. Their dropdown
values come from `docs-src/tools/helpers/schema.js`. That file is
**generated**, so do not edit it. After changing valid assays, modalities,
chemistries or index kits, regenerate it:

```bash
python3 assets/generate_tool_schema.py
```

That reads `lib/samplesheet.nf`, `lib/chemistry.nf` and `assets/indexes/`, so
the tools cannot offer something the pipeline rejects.

Navigation links are in `NAV_LINKS` in `helpers/functions.js`. They are
relative paths into the built site, so renaming a page means editing them.

## Tests and docs

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
bash tests/run_all.sh
```

The first parses every file. The second runs the self-checks: index kits,
chemistry lookup, samplesheet parsing, OverrideCycles generation, FASTQ staging
and routing, demux QC, the cellranger MultiQC table, publish layout, and a
`main.nf` compile check. Neither needs a cluster.

GitHub Actions builds this site on every push to `main` touching `docs-src/` or
`mkdocs.yml`. Pull requests build without deploying, and `mkdocs build
--strict` fails on a broken link or a page missing from the nav. Preview
locally:

```bash
pip install -r docs-src/requirements.txt
mkdocs serve
```
