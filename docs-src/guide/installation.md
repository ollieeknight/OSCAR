# Installation

## Requirements

You need Nextflow 24.04 or newer, Apptainer on every compute node, and a SLURM
scheduler. OSCAR pulls its own container images, so you do not install
cellranger, kallisto, or the rest yourself.

!!! note "Two containers you build yourself"

    OSCAR pulls everything else, but `cyto` and `cellbender` are local `.sif`
    files. Build them once and point the config at them.

```bash
apptainer build /path/to/cache/cyto.sif OSCAR/containers/cyto.def
```

Set `container_cyto` and `container_cellbender` to those paths in
`nextflow.config`, or override them on the command line.

## Get the pipeline

```bash
git clone https://github.com/ollieeknight/OSCAR
cd OSCAR
```

## Build references

`assets/references/` holds build scripts for human and mouse:

```bash
bash assets/references/human.sh          # GRCh38, cellranger + ARC
bash assets/references/mouse.sh          # GRCm38
bash assets/references/human_flex.sh     # Flex probe set
bash assets/references/human_kallisto.sh # ASAP ADT counting
```

Each script writes to the paths in `nextflow.config`. Change the `ref_*`
parameters if you put them elsewhere.

## Point the config at your cluster

`nextflow.config` ships with absolute paths for the Romagnani Lab cluster.
Change these before your first run:

| Setting | What it is |
|---------|-----------|
| `workDir` | Nextflow scratch space, needs to be fast and large |
| `apptainer.cacheDir` | Where pulled images live |
| `apptainer.runOptions` | Bind mounts for every filesystem you read or write |
| `process.queue` | Your SLURM partition, and the `gpu` queue for cellbender |
| `ref_*`, `snp_vcf`, `atac_whitelist` | Reference locations |

Cellbender needs a GPU. The `slurm` profile requests one through
`clusterOptions = '--gres=gpu:1'` on the `gpu` queue. Change both if your
scheduler names them differently.

## Check the install

Run the library self-checks. They need no cluster and no data:

```bash
nextflow run tests/test_lib.nf
```

You should see `OK: all lib/ self-checks passed`. That confirms index kit
loading, chemistry lookup, and cellranger config generation all work.
