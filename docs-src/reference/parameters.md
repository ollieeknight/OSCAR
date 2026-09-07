# Parameters

Use `--name value`. Always use underscores: `--run_until`, not `--run-until`.

## Most-used parameters

| Parameter | Use |
|---|---|
| `samplesheet` | Metadata CSV. Required. |
| `bcl_dir` | BCL folder. Required for new run. |
| `outdir` | Results folder. Default: `results`. |
| `run_name` | Results prefix. Default: BCL folder name. |
| `adt_files_dir` | ADT/HTO reference folder when not beside metadata. |
| `run_until` | Stop after `FASTQ` or `cellranger`. |
| `from_fastq`, `fastq_dir` | Start counting from existing FASTQs. |
| `from_cellranger`, `outs_dir` | Run QC on existing Cell Ranger output. |

`from_fastq` and `from_cellranger` cannot be combined.

## Several flowcells

| Parameter | Use |
|---|---|
| `extra_bcl_dirs` | Comma-separated additional BCL folders. |
| `extra_samplesheets` | Matching metadata CSVs, same order. |

Omit `extra_samplesheets` when all flowcells share primary metadata CSV.

## Optional analyses

| Parameter | Use |
|---|---|
| `run_velocity` | Set `true` for RNA velocity. Not available for Flex. |
| `viral_piscem_index` | Enables viral transcript detection for human libraries. |
| `flex_backend` | `cellranger`, `cyto`, or `both`; default `cellranger`. |
| `flex_probe_set_custom` | Semicolon-delimited custom Flex probes. |

Cluster administrators configure references, containers, and resources in `nextflow.config`.
