# OSCAR

Ollie's Single Cell Analysis for the Romagnani Lab.

OSCAR takes a raw BCL folder and a samplesheet, and gives you demultiplexed
FASTQs, count matrices, and QC output. It runs on SLURM through Nextflow DSL2,
with every tool in an Apptainer container.

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

## What runs

BCL Convert demultiplexes each lane, then falco and MultiQC report on read
quality. Counting depends on the assay. QC then runs cellbender for ambient RNA
removal and scrublet for doublet detection on GEX, AMULET and mgatk2 and MACS3
on ATAC. When a library has more than one donor, cellsnp-lite and vireo assign
cells to donors.

Two optional steps stay off by default: viral transcript detection, and RNA
velocity quantification with simpleaf.

## Start here

- [Installation](guide/installation.md) covers what you need before the first run.
- [Quickstart](guide/quickstart.md) walks through one BCL folder end to end.
- [Samplesheet](guide/samplesheet.md) documents every column and its accepted values.
- [Parameters](reference/parameters.md) lists every parameter with its default.
- [Troubleshooting](guide/troubleshooting.md) covers the failures you are most likely to hit.

Questions go to `oliver.knight@charite.de`.
