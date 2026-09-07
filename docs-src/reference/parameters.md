# Parameters

Set these on the command line with `--name value`, or edit `nextflow.config`.

!!! danger "Use underscores, not hyphens"

    Nextflow converts `--kebab-case` to `camelCase`, not to `snake_case`.
    Writing `--run-until` sets `params.runUntil`, which OSCAR never reads, so
    your run proceeds as though you passed nothing. No error, no warning.

    Use the underscore spellings below.

## Input

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samplesheet` | none | Path to metadata CSV. Required. |
| `bcl_dir` | none | BCL folder holding `RunInfo.xml`. Required unless starting downstream. |
| `outdir` | `results` | Output root. |
| `run_name` | none | Output prefix. Derived from `bcl_dir` when unset. |
| `adt_files_dir` | none | Fallback directory for ADT feature CSVs. |

## Entry points

| Parameter | Default | Description |
|-----------|---------|-------------|
| `from_fastq` | `false` | Skip demultiplexing. Needs `fastq_dir`. |
| `fastq_dir` | none | Existing FASTQs, matched by filename prefix. |
| `from_cellranger` | `false` | Run QC only. Needs `outs_dir`. |
| `outs_dir` | none | Existing `{run_name}_outs` directory. |
| `run_until` | none | `FASTQ` or `cellranger`. |

`from_fastq` and `from_cellranger` conflict. `run_until` does nothing alongside
`from_cellranger`.

## Multiple flowcells

| Parameter | Default | Description |
|-----------|---------|-------------|
| `extra_bcl_dirs` | none | Comma-separated additional BCL folders. |
| `extra_samplesheets` | none | Comma-separated samplesheets, one per extra folder. |

Supply both and the lists must match in length and order. Supply neither and
OSCAR reuses the primary samplesheet for every flowcell.

## References

| Parameter | Description |
|-----------|-------------|
| `ref_human`, `ref_mouse` | cellranger transcriptome, ARC-compatible |
| `ref_vdj_human`, `ref_vdj_mouse` | IMGT VDJ reference |
| `snp_vcf` | 1000 Genomes VCF for cellsnp-lite |
| `atac_whitelist` | ATAC barcode whitelist for ASAP ADT counting |
| `tenx_barcodes_dir` | cellranger barcode directory |
| `sequencer` | `novaseq_x` or `novaseq6000`, used when `RunInfo.xml` is unreadable |

## Flex

| Parameter | Default | Description |
|-----------|---------|-------------|
| `flex_backend` | `cellranger` | `cellranger`, `cyto`, or `both`. |
| `flex_probe_set` | 10x human v2 probe CSV | Standard probe set. |
| `flex_probe_set_custom` | none | Custom probes, semicolon-delimited. |
| `flex_samples_file` | none | Sample-to-probe-barcode CSV. Required for multiplexed Flex. |
| `flex_cyto_memory_limit` | `8GiB` | Passed to `cyto --memory-limit`. |

Give `flex_probe_set_custom` a file and OSCAR merges it with the standard set,
then hands the merged CSV to both backends. Custom probes override standard
ones by `probe_id`. The custom file looks like:

```csv
gene_id;probe_seq;probe_id;included;region;gene_name
```

## Optional analyses

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_velocity` | `false` | simpleaf spliced/unspliced quantification. |
| `spliceu_index_human`, `spliceu_index_mouse` | cluster paths | simpleaf spliceu indexes. |
| `viral_piscem_index` | cluster path | piscem index. Setting it enables viral detection. |
| `viral_t2g` | cluster path | Transcript-to-gene map for the viral index. |
| `bamtofastq_bin` | cluster path | 10x `bamtofastq` binary. |
| `macs3_qvalue` | `0.05` | MACS3 FDR threshold. |

Velocity skips Flex libraries. Viral detection runs on human libraries only.

## Containers

Every process names its image through a `container_*` parameter. Most point at
public registries and pull on first use. Two are local `.sif` files you build
yourself:

| Parameter | Default |
|-----------|---------|
| `container_cyto` | `{cache}/cyto.sif` |
| `container_cellbender` | `{cache}/cellbender_sc1.sif` |

The rest cover bclconvert, cellranger, cellranger-atac, amulet, mgatk, asap,
scrublet, falco, pigz, multiqc, kallisto, bustools, cellsnp, vireo, simpleaf,
macs3, and python. Read `nextflow.config` for the pinned versions.

!!! danger "Apptainer must be enabled"

    The `slurm` profile sets `apptainer.enabled = true`. Without it Nextflow
    ignores every `container` directive and runs bare commands on the node,
    which fails with `command not found`.

## Resources

The `slurm` profile sets cpus, memory, and walltime per process. Two need
cluster-specific attention:

- `CELLBENDER` requests the `gpu` queue with `--gres=gpu:1`
- `VIRAL_DETECT` and `SIMPLEAF_VELOCITY` pin AMD Genoa or Turin nodes

`BCLCONVERT` retries twice with doubled resources. Everything else uses
`errorStrategy = 'finish'`, so a failure lets running jobs drain and stops the
pipeline.
