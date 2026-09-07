# Architecture

For pipeline maintainers. Lab users: start with [Quickstart](../guide/quickstart.md).

## What OSCAR does

1. Validates metadata and expands sequencing indexes.
2. Demultiplexes BCL data and creates read QC reports.
3. Counts each library with assay-appropriate tool.
4. Runs applicable GEX, ATAC, and donor-assignment QC.

## Repository map

| Path | Purpose |
|---|---|
| `main.nf` | Workflow routing |
| `nextflow.config` | Profiles, parameters, resources |
| `lib/` | Metadata, chemistry, index, config helpers |
| `modules/` | Individual tool processes |
| `subworkflows/` | Process chains |
| `assets/` | Indexes, references, examples |
| `tests/test_lib.nf` | Library self-checks |

Chemistry definitions: `lib/chemistry.nf`. Read masks: `modules/demux.nf`.

```bash
nextflow run tests/test_lib.nf
```
