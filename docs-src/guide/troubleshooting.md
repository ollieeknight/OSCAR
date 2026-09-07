# Troubleshooting

## `bcl-convert: command not found`, exit 127

Nextflow ran the command outside its container. Check that your profile sets:

```groovy
apptainer.enabled = true
```

Without it Nextflow ignores every `container` directive and runs bare commands
on the compute node. Any tool can fail this way, so the same fix applies to
`cellranger: command not found` and the rest.

`BCLCONVERT` retries twice, so you see this error 3 times per lane before the
pipeline stops.

## A flag you passed did nothing

Nextflow maps `--kebab-case` to `camelCase`. Passing `--run-until FASTQ` sets
`params.runUntil`, which OSCAR never reads, so your run continues past the stage
you meant to stop at. Use `--run_until`.

The same trap catches `--from-fastq` and `--from-cellranger`.

## `unrecognised chemistry`

OSCAR matches the `chemistry` column against `lib/chemistry.nf` exactly. Both
`ARC-v1` and `ARCv1` work, as do `SC5P` and `SC5P-R2`, because both spellings
are registered. Anything else fails.

Add a chemistry by adding one entry to `chemistry_registry()`. Every consumer
reads from that map, so one edit covers demultiplexing, viral detection,
velocity, and the Flex barcode files.

## `Cannot determine OverrideCycles`

Your assay, chemistry, index type, and modality combination has no read mask.
The error names all four. Masks live in `get_override_cycles()` in
`modules/demux.nf`.

Hitting this on a valid combination means the mask table is missing an entry
rather than your samplesheet being wrong.

## `No lanes with cbcl data found`

OSCAR looks for `{bcl_dir}/Data/Intensities/BaseCalls/L00*/C1.1/*.cbcl`. An
incomplete transfer off the sequencer is the usual cause. Confirm the run
finished copying before you start.

## `Demux group has mixed index lengths`

One demultiplexing group holds both 8bp and 10bp indexes, and BCL Convert
cannot mask both in a single sample sheet. Groups key on assay, index type,
chemistry, modality, and BCL folder, so this means one library mixes index
lengths within a modality. Correct the `index` column.

## ADT libraries fail at `cellranger multi`

The `[feature]` reference never resolved. OSCAR warns at parse time, then fails
at counting. Check that `adt_file` names a real CSV, then that the CSV sits in
one of the three search locations from
[Samplesheet](samplesheet.md#adt-and-hto-feature-references).

## Custom Flex probes missing from output

Earlier releases generated the merged probe CSV and then handed cellranger the
standard one, so custom probes reached the cyto backend alone. Current versions
pass the merged CSV to both. Confirm your `--flex_probe_set_custom` path exists, then
check `multi_config.csv` in the work directory for a `probe-set,` line that
points at `merged_probes_cellranger.csv`.

## Velocity skipped a library

Velocity needs a chemistry with a `velocity` entry in the registry. Flex is
probe-based and carries no intronic signal, so OSCAR filters those libraries out
without failing. ATAC and `NA` behave the same way.

## Reading a failure

Nextflow prints the work directory for the failed task. Go there:

```bash
cd /path/to/work/3e/b0a279...
cat .command.sh    # the exact command
cat .command.err   # stderr
cat .command.log   # both streams
```

!!! tip

    `.command.sh` shows the command after Nextflow substituted every variable,
    so it tells you what ran rather than what the module says.

Resume from the last successful step after fixing the cause:

```bash
nextflow run main.nf -resume ...
```

## Checking pipeline logic without a cluster

```bash
nextflow lint main.nf lib modules subworkflows nextflow.config
nextflow run tests/test_lib.nf
```

The first parses every file. The second exercises index loading, chemistry
lookup, samplesheet parsing, and cellranger config generation against the
shipped example samplesheet. Neither needs SLURM or data.
