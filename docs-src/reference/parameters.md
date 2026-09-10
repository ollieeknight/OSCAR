# Parameters

Use `--name value`. Always use underscores: `--run_from`, not `--run-from`.

## Most-used parameters

| Parameter | Use |
|---|---|
| `samplesheet` | Metadata CSV. Required. |
| `bcl_dir` | BCL folder. Required for new run. |
| `outdir` | Results folder. Default: `results`. |
| `run_name` | Results prefix. Default: BCL folder name. |
| `adt_files_dir` | ADT/HTO reference folder when not beside metadata. |
| `run_from` | Entry point: `bcl` (default), `fastq`, or `cellranger`. |
| `fastq_dir` | Existing FASTQs. Required for `--run_from fastq`. |
| `outs_dir` | Existing Cell Ranger output. Required for `--run_from cellranger`. |
| `extras` | Optional analyses, comma-separated: `velocity`, `viral`. |

## Several flowcells

| Parameter | Use |
|---|---|
| `extra_bcl_dirs` | Comma-separated additional BCL folders. |
| `extra_samplesheets` | Matching metadata CSVs, same order. |

Omit `extra_samplesheets` when all flowcells share primary metadata CSV.

## Optional analyses

| Parameter | Use |
|---|---|
| `extras` | `velocity` for RNA velocity (not available for Flex), `viral` for viral transcript detection (human libraries). Combine as `--extras velocity,viral`. |
| `viral_piscem_index`, `viral_t2g`, `bamtofastq_bin` | Reference paths for `--extras viral`. |
| `spliceu_index_human`, `spliceu_index_mouse` | Reference paths for `--extras velocity`. |
| `flex_backend` | `cellranger`, `cyto`, or `both`; default `cellranger`. |
| `flex_probe_set_custom` | Semicolon-delimited custom Flex probes. |

Cluster administrators configure references, containers, and resources in `nextflow.config`.
