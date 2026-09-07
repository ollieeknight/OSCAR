# Samplesheet

One CSV, one row per library-modality pair. A DOGMA library sequenced for GEX,
ATAC, and ADT takes three rows.

Build one with the [metadata generator](../tools/metadata_generator.html), or
copy `assets/example_metadata.csv`.

## Columns

All eleven are required. Extra columns are ignored.

| Column | Accepts | Notes |
|--------|---------|-------|
| `assay` | GEX, CITE, DOGMA, ATAC, Multiome, ASAP, Flex | Case-insensitive |
| `experiment_id` | free text | No commas |
| `historical_number` | integer | Your running experiment count |
| `replicate` | free text | Usually A, B, C |
| `modality` | GEX, ATAC, ADT, HTO, VDJ-T, VDJ-B, CRISPR, GENO | Case-sensitive |
| `chemistry` | see below | Case-sensitive |
| `index_type` | SI, DI, NA | |
| `index` | kit code or raw sequence | `SI-TT-A1`, or `CATAGCCG` |
| `species` | human, mouse | Case-insensitive |
| `n_donors` | integer or NA | Above 1 triggers vireo |
| `adt_file` | basename or NA | Resolves to `{adt_file}.csv` |

## Chemistry values

| Value | Kit |
|-------|-----|
| `SC3Pv2` | 3' v2 |
| `SC3Pv3` | 3' v3 |
| `SC3Pv4` | 3' v4 (GEM-X) |
| `SC5P`, `SC5P-R2` | 5' v2 |
| `SC5Pv3` | 5' v3 |
| `ARCv1`, `ARC-v1` | Multiome and DOGMA |
| `Flex-v2-R1`, `Flex-v2-RNA-R2` | Fixed RNA Profiling v2 |
| `ATAC` | ATAC and ASAP |
| `NA` | Modalities without their own chemistry |

Both spellings of the Multiome and 5' v2 chemistries work. OSCAR validates this
column against `lib/chemistry.nf`, which is where you add a new chemistry.

## How rows become library IDs

OSCAR builds two identifiers from every row:

```
library_id = {assay}_{experiment_id}_exp{historical_number}_lib{replicate}
id         = {library_id}_{modality}
```

The row `DOGMA,Tissue_HSPC,1,A,GEX,ARC-v1,DI,SI-TT-A1,human,2,Tissue_HSPC_ADT`
produces `DOGMA_Tissue_HSPC_exp1_libA` and
`DOGMA_Tissue_HSPC_exp1_libA_GEX`.

Rows sharing a `library_id` get counted together, which is how OSCAR knows that
your GEX, ADT, and HTO rows belong in one `cellranger multi` run. Output folders
use `library_id`, so keep those four columns consistent across a library's rows.

## Indexes

Name a 10x kit and OSCAR expands it from `assets/indexes/`:

- `SI-GA-*`, `SI-NA-*` are single-index and expand to four sequences
- `SI-TT-*`, `SI-TN-*`, `SI-TS-*` are dual-index

Anything unrecognised is treated as a raw i7 sequence, which is how the 8bp
TruSeq indexes on ADT and HTO libraries work (`CATAGCCG` above).

OSCAR picks the i5 orientation for you by reading `RunInfo.xml`. Instrument
`VH` and `NDX` and `LH` and `FS` get forward i5; `A` and `NB` and `NS` and `MN`
get reverse complement. Without a readable `RunInfo.xml`, OSCAR falls back to
`--sequencer` and warns.

## Spacer rows

Comma-only rows separate experiments visually and OSCAR skips them:

```csv
DOGMA,Tissue_HSPC,1,A,ADT,ARC-v1,DI,CATAGCCG,human,2,Tissue_HSPC_ADT
,,,,,,,,,,
CITE,Helicobacter_colitis,2,A,GEX,SC3Pv3,SI,SI-GA-A1,mouse,NA,Hh_colitis
```

## ADT and HTO feature references

`adt_file` names a CSV without its extension. OSCAR looks in three places and
takes the first hit:

1. `{samplesheet_dir}/adt_files/{adt_file}.csv`
2. `{samplesheet_dir}/../adt_files/{adt_file}.csv`
3. `{adt_files_dir}/{adt_file}.csv`

Find none and OSCAR warns, then `cellranger multi` fails on the missing
`[feature]` reference. The file follows the 10x feature reference format:

```csv
id,name,read,pattern,sequence,feature_type
CD3,CD3,R2,5PNNNNNNNNNN(BC),CTCATTGTAACTCCT,Antibody Capture
```

`assets/example_adt_totalseq_{a,b,c,d}.csv` cover the TotalSeq panels, or build
one with the [feature barcode generator](../tools/adt_generator.html).

## Donor demultiplexing

Set `n_donors` above 1 and OSCAR runs cellsnp-lite then vireo for that library.
This needs human data, since `snp_vcf` points at a 1000 Genomes VCF. Mouse
libraries skip it whatever you put in the column.

Use `NA` or `1` for single-donor libraries.
