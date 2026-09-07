# OSCAR

OSCAR processes a sequencing run into FASTQs, count matrices, and quality-control results.

## Prepare these before you run

1. **Sequencing data:** a completed BCL folder.
2. **Metadata:** make one CSV for the run with the [metadata generator](tools/metadata_generator.html).
3. **Feature barcodes:** for ADT or HTO libraries, make a reference CSV with the [feature barcode generator](tools/adt_generator.html). Put it in `adt_files/` beside the metadata CSV.

Then copy the command from [Quickstart](guide/quickstart.md).

## Supported assays

| Assay | Library types |
|---|---|
| GEX | GEX, VDJ-T, VDJ-B, CRISPR |
| CITE | GEX, ADT, HTO |
| Flex | GEX |
| ATAC | ATAC |
| Multiome | GEX, ATAC |
| DOGMA | GEX, ATAC, ADT, HTO |
| ASAP | ATAC, ADT, HTO |

## Help

- [Quickstart](guide/quickstart.md): run one flowcell.
- [Samplesheet](guide/samplesheet.md): metadata fields and feature references.
- [Installation](guide/installation.md): set up OSCAR once.
- [Troubleshooting](guide/troubleshooting.md): fix common failures.

Questions: `oliver.knight@charite.de`.
