# Architecture

## Layout

```
main.nf              routing and channel logic
nextflow.config      params, profiles, per-process resources
lib/                 pure functions, no processes
├── chemistry.nf     chemistry registry
├── indexes.nf       index kits, sequencer detection
├── samplesheet.nf   parsing and validation
└── multi_config.nf  cellranger multi config generation
modules/             one process per tool
subworkflows/        process chains
tests/test_lib.nf    self-checks for lib/
assets/
├── indexes/         10x index kit CSVs
├── references/      reference build scripts
└── example_*.csv    samplesheet and ADT templates
```

Nothing in `lib/` declares a process or reads `params`. You call these functions
from `tests/test_lib.nf` without a scheduler, which is what makes the self-checks
possible.

## Flow

```
samplesheet ──> parse ──> per-flowcell index resolution
                            │
                            v
                  GENERATE_SAMPLESHEET ──> BCLCONVERT (per lane)
                            │
                            v
                  VALIDATE_FASTQ, FALCO ──> MULTIQC
                            │
        ┌───────────────────┼───────────────────┐
        v                   v                   v
  CELLRANGER_MULTI   CELLRANGER_ATAC      COUNT_ADT
   GEX CITE Flex       ATAC assays       ASAP ADT/HTO
   Multiome DOGMA
        │                   │
        v                   v
     QC_GEX              QC_ATAC
  cellbender          AMULET mgatk2
  scrublet            MACS3
        │                   │
        └────── vireo ──────┘
             n_donors > 1
```

## The chemistry registry

`lib/chemistry.nf` holds one map keyed by chemistry name. Each entry carries the
demultiplexing family, the 10x barcode whitelist, the simpleaf chemistry string,
the velocity string, and the Flex barcode files.

Every consumer reads from it. Samplesheet validation derives its accepted values
from the keys, `modules/demux.nf` resolves its `OverrideCycles` key through
`get_chemistry_family()`, and the Flex modules call `get_flex_barcode_file()`
rather than matching on the chemistry name themselves.

Adding a chemistry means adding one entry. Before this consolidation the same
knowledge sat in four files, and `Flex-v2-RNA-R2` passed validation while having
no entry in any lookup table.

Accessors fail loudly on a chemistry that exists but lacks the field you asked
for, so a half-added entry breaks at the point of use instead of doing nothing.
`get_velocity_chemistry()` returns null rather than failing, because callers
filter on it.

## cellranger multi config

`multi_config.csv` gets built in two places. `lib/multi_config.nf` writes the
`[gene-expression]`, `[vdj]`, and `[feature]` sections before staging.
`modules/count_gex.nf` writes `[libraries]` in Python inside the process, since
only the process knows which staged FASTQ runs cleared the read-count filter and
what lane numbers they ended up with.

`COUNT_GEX` calls the header builder rather than `main.nf`, which is what lets
the merged custom probe CSV from `FLEX_PROBE_PREPARE` reach the `probe-set,`
line.

## Read masks

`get_override_cycles()` in `modules/demux.nf` maps assay, chemistry family,
index type, and modality to a BCL Convert `OverrideCycles` string. Masks stay
there rather than in the registry because they depend on all four, not on
chemistry alone.

`GENERATE_SAMPLESHEET` resolves the wildcards at runtime by reading cycle counts
from `RunInfo.xml`, since BCL Convert 4.x rejects `*`.

Two corrections apply. An 8bp TruSeq index turns `I10N*` into `I8N2`. A
single-index library on a 4-read dual-index flowcell clamps index 1 to the real
sequence length and masks index 2 completely.

## Multi-flowcell handling

Each BCL folder gets its own sequencer detection and its own index resolution,
so mixing instruments works. OSCAR then deduplicates library metadata by `id`
while preserving order, and groups FASTQs by `library_id` across flowcells
before counting.

Entries sort by first-file URI for deterministic flowcell ordering. Files within
a flowcell keep the name sort that BCL Convert produced, which puts R1 before R2.
A global name sort would group every R1 together across flowcells and pair the
wrong reads.
