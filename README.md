# OSCAR

**Ollie's Single Cell Analysis for the Romagnani Lab**

OSCAR turns a finished sequencing run into demultiplexed FASTQs, count
matrices, and QC reports. Nextflow DSL2 on SLURM, every tool in an Apptainer
container.

Documentation: **[ollieeknight.github.io/OSCAR](https://ollieeknight.github.io/OSCAR/)**

- Sequencing libraries? Fill in the
  [metadata generator](https://ollieeknight.github.io/OSCAR/tools/metadata_generator.html)
  and send it over. See
  [what the fields mean](https://ollieeknight.github.io/OSCAR/guide/filling-these-in/).
- Running it? [Running OSCAR](https://ollieeknight.github.io/OSCAR/guide/running/).
- Taking it over? [Maintaining OSCAR](https://ollieeknight.github.io/OSCAR/reference/maintaining/).

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

Beyond counting: fastp and MultiQC on reads, cellbender and scrublet on GEX,
AMULET, mgatk2 and MACS3 on ATAC, cellsnp-lite and vireo for multi-donor
libraries. Viral detection and RNA velocity are off by default; turn them on
with `--extras viral,velocity`.

## Running

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir     /path/to/R463_bcl
```

OSCAR is configured for the Romagnani lab's Charité SC cluster. The reference,
scratch and bind paths in `nextflow.config` resolve nowhere else; see
[Maintaining OSCAR](https://ollieeknight.github.io/OSCAR/reference/maintaining/)
to point them elsewhere.

## Checking your setup

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
bash tests/run_all.sh
```

The first parses every file, the second runs every self-check. Neither needs a
cluster, containers, or sequencing data.

## Citation and licence

MIT, see [LICENSE](LICENSE). Cite with the metadata in
[CITATION.cff](CITATION.cff).

Oliver Knight, `oliver.knight@charite.de`
