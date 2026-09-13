# Filling these in

Two files describe your run. Everyone sends a metadata CSV. If you have
antibody or hashtag libraries, you also send one feature barcode CSV per panel.
Build both here:

1. [Metadata generator](../tools/metadata_generator.html)
2. [Feature barcode generator](../tools/adt_generator.html)

Both write the file for you. The menus only offer values OSCAR accepts, so you
cannot mistype a chemistry or an index code. What the menus cannot decide for
you is how many rows you need, and that is what this page is about.

## One row per library type

Each row describes **one library type from one library**, not one sample.

A plain gene-expression library is one row. A CITE library sequenced for gene
expression and surface antibodies is two rows, one `GEX` and one `ADT`. A
DOGMA library is four.

Rows belonging to the same library must agree on these four fields:

- `assay`
- `experiment_id`
- `historical_number`
- `replicate`

OSCAR groups rows by those four. If they match, the rows are one library; if
any one of them differs, you have accidentally described two libraries.

## What each field means

| Field | What to enter |
|---|---|
| `assay` | The experiment type: GEX, CITE, ATAC, Multiome, DOGMA, ASAP, or Flex |
| `experiment_id` | A short name for the experiment. No commas |
| `historical_number` | Your experiment number |
| `replicate` | Which replicate this library is, usually A, B or C |
| `modality` | The library type this row describes |
| `chemistry` | The 10x kit used. Take it from the kit box |
| `index_type` | `SI` for a single-index kit, `DI` for dual, `NA` for none |
| `index` | The 10x index code, or the raw sequence for a custom index |
| `species` | `human` or `mouse` |
| `n_donors` | How many people the sample came from. `1` or `NA` for one |
| `adt_file` | The feature barcode filename, without `.csv`. `NA` if unused |

`n_donors` above one turns on donor demultiplexing: OSCAR separates the cells
genetically and tells you which donor each came from. Get this right. Nobody
can work it out afterwards.

## Antibody and hashtag libraries

For any `ADT` or `HTO` row, make a reference with the
[feature barcode generator](../tools/adt_generator.html). Pick the species,
the TotalSeq panel, and the markers you stained with.

Name the file after the panel, not after a sample: several libraries stained
with the same panel share one reference. Put the name in `adt_file` **without
the `.csv`**. A file called `totalseq_a_myeloid.csv` is entered as
`totalseq_a_myeloid`.

Send the reference CSVs along with the metadata CSV.

## Before you send it

- One metadata CSV per sequencing run.
- Every library type you sequenced has a row.
- Rows from one library agree on assay, experiment ID, historical number and replicate.
- Every `adt_file` name matches a reference CSV you are sending.
- `n_donors` is right for pooled samples.

Send both to `oliver.knight@charite.de`. You will get back count matrices and a
QC report.
