# Troubleshooting

## Start here

Open failed task folder reported by Nextflow:

```bash
cd /path/to/work/3e/b0a279...
cat .command.err
cat .command.sh
```

Fix cause, then rerun with `-resume`. Nextflow reuses completed work.

```bash
nextflow run main.nf -profile slurm -resume [your usual options]
```

## Common problems

### `command not found`

Apptainer did not start. Check selected profile sets:

```groovy
apptainer.enabled = true
```

### A parameter did nothing

Use underscores: `--run_until`, `--from_fastq`, and `--from_cellranger`. Hyphenated forms do not set OSCAR parameters.

### `unrecognised chemistry`

Choose chemistry from metadata generator. Values must match exactly.

### ADT or HTO fails during counting

Check `adt_file` matches real `{adt_file}.csv`. Check location in [Metadata and feature barcodes](samplesheet.md#adt-and-hto-references).

### `No lanes with cbcl data found`

BCL transfer is incomplete. Confirm sequencer run finished copying, then retry.

### `Demux group has mixed index lengths`

One library type combines 8 bp and 10 bp indexes. Correct `index` values in metadata CSV.

### Velocity missing

Velocity does not run for Flex libraries. Expected.
