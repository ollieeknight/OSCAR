# Plan 01: Remove TOP_UNKNOWN_BARCODES and publish all BCL Convert FASTQs

## Objective

1. Remove the `TOP_UNKNOWN_BARCODES` process entirely (process definition + subworkflow call + import).
2. Publish **all** FASTQs output by BCLCONVERT (named samples + Undetermined) directly from the BCLCONVERT process, replacing the current per-sample publish in VALIDATE_FASTQ.

## Context

- `modules/demux.nf` — contains process definitions: `BCLCONVERT`, `TOP_UNKNOWN_BARCODES`, `VALIDATE_FASTQ`
- `subworkflows/demux.nf` — orchestrates demux; imports and calls `TOP_UNKNOWN_BARCODES`; aggregates undetermined R1s to feed it
- Currently named-sample FASTQs are published via `VALIDATE_FASTQ.publishDir` → `${outdir}/${run_name}_fastq/`
- Undetermined FASTQs are never published; only the unknown-barcode summary CSV is
- `BCLCONVERT` runs per-lane; all outputs land in `fastqs/*.fastq.gz` inside the work dir

---

## Phase 0: Verify current publishing paths (read-only check)

Before editing, confirm these exact strings in the files:

```
grep -n "publishDir" OSCAR/modules/demux.nf
grep -n "TOP_UNKNOWN_BARCODES\|undetermined" OSCAR/subworkflows/demux.nf
```

Expected:
- `VALIDATE_FASTQ` has `publishDir { "${params.outdir}/${meta.run_name}_fastq" }, mode: 'copy'`
- `TOP_UNKNOWN_BARCODES` has `publishDir { ... "${params.outdir}/${run}_fastq" }`
- `subworkflows/demux.nf` includes `TOP_UNKNOWN_BARCODES` import and a `BCLCONVERT.out.undetermined` pipeline block

---

## Phase 1: Add publishDir to BCLCONVERT and remove publishDir from VALIDATE_FASTQ

**File:** `modules/demux.nf`

### 1a. Add publishDir to BCLCONVERT

Insert after the `tag` and `container` lines of the BCLCONVERT process:

```groovy
publishDir {
    def run = bcl_dir.name.replaceAll(/_bcl.*$/, '')
    "${params.outdir}/${run}_fastq"
}, mode: 'copy', pattern: "fastqs/*.fastq.gz", saveAs: { fn -> file(fn).name }
```

This mirrors the pattern used in `TOP_UNKNOWN_BARCODES` for deriving `run` from `bcl_dir.name`.

### 1b. Remove publishDir from VALIDATE_FASTQ

Delete the `publishDir` directive from `VALIDATE_FASTQ` — BCLCONVERT now publishes the FASTQs.
Keep the `pigz -t` validation script intact; VALIDATE_FASTQ still runs for integrity checking, just no longer publishes.

**Verification:**
```
grep -n "publishDir" OSCAR/modules/demux.nf
# Expected: only BCLCONVERT and TOP_UNKNOWN_BARCODES (latter removed in Phase 2)
```

---

## Phase 2: Remove TOP_UNKNOWN_BARCODES

### 2a. `modules/demux.nf` — delete the process block

Delete the entire `TOP_UNKNOWN_BARCODES` process (from `// ─── TOP_UNKNOWN_BARCODES` comment to closing `}`).

### 2b. `subworkflows/demux.nf` — remove import and call

Remove:
```groovy
include { TOP_UNKNOWN_BARCODES  } from '../modules/demux'
```

Remove the entire block that aggregates undetermined and calls TOP_UNKNOWN_BARCODES:
```groovy
// Aggregate undetermined R1s across lanes per demux group, count top unknown barcodes.
// Runs independently — no downstream process waits on it.
BCLCONVERT.out.undetermined
    .groupTuple(by: [0, 2])
    .map { demux_key, metas_per_lane, bcl_dir_name, r1_lists ->
        [demux_key, metas_per_lane[0], bcl_dir_name, r1_lists.flatten()]
    }
    | TOP_UNKNOWN_BARCODES
```

### 2c. `modules/demux.nf` BCLCONVERT output — remove undetermined emit

Remove from BCLCONVERT `output:`:
```groovy
tuple val(demux_key), val(metas), val(bcl_dir.name), path("fastqs/Undetermined_S*_L*_R1_*.fastq.gz"), optional: true, emit: undetermined
```

The `fastqs/*.fastq.gz` pattern in the new publishDir already covers Undetermined files.

**Verification:**
```
grep -rn "TOP_UNKNOWN_BARCODES\|undetermined" OSCAR/modules/demux.nf OSCAR/subworkflows/demux.nf
# Expected: zero matches
```

---

## Phase 3: Final verification

```bash
# No references to removed process remain
grep -rn "TOP_UNKNOWN_BARCODES" OSCAR/

# BCLCONVERT publishDir present
grep -A4 "process BCLCONVERT" OSCAR/modules/demux.nf | grep publishDir

# VALIDATE_FASTQ publishDir absent
grep -A10 "process VALIDATE_FASTQ" OSCAR/modules/demux.nf | grep -v publishDir

# Pipeline still parses (requires nextflow)
nextflow -C OSCAR/nextflow.config run OSCAR/main.nf --help 2>&1 | head -5
```

---

## Summary of changes

| File | Change |
|------|--------|
| `modules/demux.nf` | Add `publishDir` to BCLCONVERT; remove `publishDir` from VALIDATE_FASTQ; delete `TOP_UNKNOWN_BARCODES` process; remove `undetermined` emit from BCLCONVERT |
| `subworkflows/demux.nf` | Remove `TOP_UNKNOWN_BARCODES` import; remove undetermined aggregation + call block |
