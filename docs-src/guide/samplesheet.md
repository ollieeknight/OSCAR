# Metadata and feature barcodes

Create the metadata CSV before running OSCAR. Use the [metadata generator](../tools/metadata_generator.html); it writes required headings and accepted values.

## One row per library type

Each row is one modality from one library. A CITE library with GEX and ADT needs two rows. Keep `assay`, `experiment_id`, `historical_number`, and `replicate` identical for rows from the same library.

## Required fields

| Field | Enter |
|---|---|
| `assay` | GEX, CITE, DOGMA, ATAC, Multiome, ASAP, or Flex |
| `experiment_id` | Short experiment name; no commas |
| `historical_number` | Experiment number |
| `replicate` | Library replicate, usually A, B, or C |
| `modality` | GEX, ATAC, ADT, HTO, VDJ-T, VDJ-B, CRISPR, or GENO |
| `chemistry` | Chemistry shown in generator |
| `index_type` | `SI`, `DI`, or `NA` |
| `index` | 10x kit code, or raw index sequence |
| `species` | `human` or `mouse` |
| `n_donors` | Number of human donors; `NA` or `1` for one donor |
| `adt_file` | Reference filename without `.csv`; `NA` when unused |

## ADT and HTO references

For ADT or HTO, use the [feature barcode generator](../tools/adt_generator.html) to select species, TotalSeq panel, and output format.

Save `{adt_file}.csv` in one of these locations:

1. `adt_files/` beside `metadata.csv`.
2. `adt_files/` in metadata folder's parent directory.
3. Directory given with `--adt_files_dir`.

`adt_file` must match filename exactly, without `.csv`.

## Check before running

- Use a separate metadata CSV for each sequencing run.
- Use exact chemistry and modality values offered by generator.
- Check every ADT/HTO filename exists.
- OSCAR ignores comma-only spacer rows.

OSCAR counts rows with matching first four fields as one library.
