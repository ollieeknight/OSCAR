# OSCAR

OSCAR turns a finished sequencing run into count matrices and quality-control
reports. It is the Romagnani lab's single-cell processing pipeline.

## Who does what

**If you have libraries to sequence,** you need two files from this site and
nothing else. Fill in the [metadata generator](tools/metadata_generator.html).
If you have antibody or hashtag libraries, also fill in the
[feature barcode generator](tools/adt_generator.html). Send both to Ollie.
[Filling these in](guide/filling-these-in.md) explains the fields.

You do not need cluster access, and you do not run OSCAR yourself.

**If you are running OSCAR,** see [Running OSCAR](guide/running.md).

**If you are taking OSCAR over,** see [Maintaining OSCAR](reference/maintaining.md).

## Two clusters

The lab uses two separate machines, and they are easy to confuse.

| | Charité SC | BIH-CUBI |
|---|---|---|
| Used for | running OSCAR | RStudio, Jupyter, your analysis |
| Set up by | Ollie | [the HPC cluster guide](https://github.com/Romagnani-Lab/bih-cubi-romagnani) |

OSCAR runs on Charité SC. Its results are copied to wherever you analyse them.
Everything in the HPC cluster guide is about BIH-CUBI and does not apply here.

## What OSCAR supports

| Assay | Library types |
|---|---|
| GEX | GEX, VDJ-T, VDJ-B, CRISPR |
| CITE | GEX, ADT, HTO |
| Flex | GEX |
| ATAC | ATAC |
| Multiome | GEX, ATAC |
| DOGMA | GEX, ATAC, ADT, HTO |
| ASAP | ATAC, ADT, HTO |

Every run also gets read-level QC. GEX libraries get ambient-RNA and doublet
calls; ATAC libraries get doublet, mitochondrial and peak calls. For any library
with more than one donor, OSCAR assigns each cell back to its donor by
genotype.

Questions: `oliver.knight@charite.de`.
