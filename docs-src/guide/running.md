# Running OSCAR

On Charité SC, with a metadata CSV and a finished BCL folder.

```bash
nextflow run main.nf -profile slurm \
  --samplesheet /path/to/metadata.csv \
  --bcl_dir     /path/to/R463_bcl
```

Results land beside the flowcell: FASTQs and read QC in `{run}_fastq`, counts
and QC in `{run}_outs`, both siblings of `{run}_bcl`. Open
`{run}_fastq/multiqc/multiqc_report.html` first.

```
{outdir}/
├── {run_name}_fastq/multiqc/multiqc_report.html
└── {run_name}_outs/
    ├── {library_id}/outs/        GEX, CITE, Flex, Multiome, DOGMA
    └── {library_id}_ATAC/outs/   ATAC, Multiome, DOGMA, ASAP
```

OSCAR looks for feature barcode CSVs in `adt_files/` beside the metadata CSV,
then in `adt_files/` one directory up. Point `--adt_files_dir` at them when they
live somewhere else.

## When a run fails

Nextflow prints the failed task's work directory. Read its error, fix the
cause, and rerun the same command with `-resume`. Nextflow reuses completed
work.

```bash
cd /path/to/work/3e/b0a279...
cat .command.err
```

Resource kills (out of memory, walltime) retry automatically with more memory
and time, twice, before the run stops.

## Options

| Option | When you need it |
|---|---|
| `--run_name` | Name the output something other than the BCL folder |
| `--outdir` | Write results somewhere other than beside the flowcell |
| `--adt_files_dir` | Feature barcode CSVs kept outside the metadata folder |
| `--run_from fastq` with `--fastq_dir` | Skip demultiplexing |
| `--run_from cellranger` with `--outs_dir` | Rerun QC on existing counts |
| `--extras velocity,viral` | RNA velocity and viral detection, both off by default |
| `--flex_backend` | Flex runs: `cellranger`, `cyto`, or `both` |

Parameter names use underscores. `--run-from` sets nothing.

`velocity` does not run for Flex libraries. `viral` is human only.

## One library across several flowcells

```bash
--extra_bcl_dirs     /path/to/R462_bcl,/path/to/R463_bcl \
--extra_samplesheets /path/to/R462.csv,/path/to/R463.csv
```

Both lists in the same order. Omit `--extra_samplesheets` when every flowcell
shares one metadata CSV. Each flowcell's FASTQs stay beside its own BCL folder;
counts collect under a single `{run}_outs` named by `--run_name`.

## Tuning

Rarely needed; defaults are set for our runs.

| Parameter | Default | Effect |
|---|---|---|
| `macs3_qvalue` | `0.05` | ATAC peak-calling FDR threshold |
| `demux_qc_dropout_ratio` | `100` | How far below its lane's median a library must fall to be reported as a dropout |
| `demux_qc_unknown_pct` | `1` | Unknown barcode at or above this percent of lane reads is reported |
| `flex_cyto_memory_limit` | `8GiB` | Memory ceiling passed to cyto |
| `sequencer` | `novaseq_x` | Fallback only. Read from `RunInfo.xml` when present |

Demux QC warnings are advisory and never fail a run.

## Checking the pipeline still works

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
bash tests/run_all.sh
```

Neither needs a cluster, containers, or sequencing data.
