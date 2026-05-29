# SPRING FASTQ Compression Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** (1) Auto-detect `.spring` files when `--from_fastq` is used and decompress them before the pipeline runs. (2) After BCL demux, compress R1+R2 FASTQ pairs to SPRING format for archival instead of copying plain FASTQ files.

**Architecture:**
Two independent features using the same container and process patterns. `SPRING_DECOMPRESS` converts `.spring` → `.fastq.gz` pairs inside the work dir before FASTQ_QC; output is mixed with regular `ch_fastqs`. `SPRING_COMPRESS` is added after `VALIDATE_FASTQ` in the FASTQ_QC subworkflow: R1/R2 pairs are grouped and compressed; VALIDATE_FASTQ's publishDir is modified to skip R1/R2 files (returning `null` from `saveAs` suppresses publishing). Index reads (I1/I2) remain published as-is.

**Tech Stack:** Nextflow DSL2, spring 1.1.1, Apptainer. Container SIF already cached at `/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-spring_1.1.1--h4ac6f70_3.img`.

---

## Phase 0: Confirmed facts

### Spring 1.1.1 CLI (from upstream README)

```bash
# Compress paired-end (single .spring output)
spring -c -i r1.fastq.gz r2.fastq.gz -o sample.spring -g -t <N>

# Decompress paired-end → two explicit .fastq.gz files
spring -d -i sample.spring -o r1.fastq.gz r2.fastq.gz -g -t <N>
```

Key facts:
- `-g` flag required when input or output is gzipped FASTQ
- Paired-end: ONE `.spring` file encodes both R1+R2
- Output naming on decompress: must give two explicit output filenames; without them, spring appends `.1`/`.2` to a single name (undesirable)
- Detection by extension: `.spring` — the file is internally a TAR archive
- Thread flag: `-t <N>`

### OSCAR FASTQ channel structure (from code)

`VALIDATE_FASTQ` (modules/demux.nf:316–333):
```groovy
input:  tuple val(meta), val(fastq_dir), path(fastq), val(fastq_name)
output: tuple val(meta), val(fastq_dir), path(fastq), emit: fastq
publishDir { "${params.outdir}/${meta.run_name}_fastq" }, mode: 'copy',
    saveAs: { fn -> file(fn).name }
```

`FASTQ_QC` subworkflow (subworkflows/fastq_qc.nf) input/output:
```
input:  ch_fastqs  = [meta, fastq_dir_str, [fq_files]]
emit:   fastqs     = [meta, [unique_fq_dirs], [validated_fqs]]
        falco_reports
```

`--from_fastq` path in main.nf (lines 341–354):
- Glob: `${params.fastq_dir}/**/${meta.id}*.fastq.gz`
- Groups by parent directory
- Emits `ch_fastqs: [meta, fastq_dir_str, [matched_fqs]]`

**Anti-patterns:**
- Do not pass spring decompression output path as a `val` string — Nextflow won't know the work dir path at channel-build time; emit `path()` outputs and derive `fq_dir` from `.parent` in a `.map{}`
- Do not use `.any {}`, `.collect {}`, `.flatten()` on groupTuple output inside a process script block — use `.each { fl << it }` materialisation in `.map{}` operators
- Do not use `saveAs: { fn -> null }` as a boolean — it must explicitly return `null` (String) for Nextflow to suppress publishing
- Do not compress I1/I2 index reads with spring — spring is for read pairs; index reads are tiny and not useful post-demux

---

## Task 1: Add `container_spring` param to nextflow.config

**Files:**
- Modify: `OSCAR/nextflow.config`

- [ ] **Step 1: Add the spring container param**

In `nextflow.config`, add after `container_pigz`:
```groovy
container_spring     = "/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-spring_1.1.1--h4ac6f70_3.img"
```

- [ ] **Step 2: Commit**

```bash
git add OSCAR/nextflow.config
git commit -m "config: add container_spring param for spring 1.1.1 FASTQ compression"
```

---

## Task 2: Add SPRING_DECOMPRESS process to modules/demux.nf

**Files:**
- Modify: `OSCAR/modules/demux.nf` — add process after the VALIDATE_FASTQ block (~line 334)

The process:
- Input: one `.spring` file per paired library
- Output: `r1.fastq.gz` + `r2.fastq.gz` in the work directory
- Naming: derive output names from the spring filename (strip `.spring` suffix, add `_R1_001.fastq.gz` / `_R2_001.fastq.gz`)

- [ ] **Step 1: Add SPRING_DECOMPRESS process**

Insert after VALIDATE_FASTQ (after line 333 of modules/demux.nf):

```groovy
// ─── SPRING_DECOMPRESS ───────────────────────────────────────────────────────
// Decompress a paired-end .spring archive → R1 + R2 .fastq.gz files.
// Called in --from_fastq mode when spring files are detected in params.fastq_dir.

process SPRING_DECOMPRESS {
    tag "$meta.id"
    label 'process_medium'
    container "${params.container_spring}"

    input:
    tuple val(meta), path(spring_file)

    output:
    tuple val(meta), path("*_R1_001.fastq.gz"), path("*_R2_001.fastq.gz"), emit: fastqs
    path "versions.yml", emit: versions

    script:
    def prefix = spring_file.name.replaceAll(/\.spring$/, '')
    """
    spring -d \\
        -i ${spring_file} \\
        -o ${prefix}_R1_001.fastq.gz ${prefix}_R2_001.fastq.gz \\
        -g \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spring: \$(spring --version 2>&1 | head -1 | sed 's/Spring v//')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Verify the process compiles (syntax check)**

```bash
# On cluster:
nextflow run OSCAR/main.nf --help
# Expected: no groovy compilation errors
```

- [ ] **Step 3: Commit**

```bash
git add OSCAR/modules/demux.nf
git commit -m "feat: add SPRING_DECOMPRESS process for --from_fastq spring file detection"
```

---

## Task 3: Wire SPRING_DECOMPRESS into main.nf `--from_fastq` path

**Files:**
- Modify: `OSCAR/main.nf` — the `--from_fastq` block (~lines 341–354) and include statement

- [ ] **Step 1: Add SPRING_DECOMPRESS to the include statement**

Find in main.nf:
```groovy
include { VIRAL_DETECT }   from './modules/qc'
```
Or wherever demux modules are included. Add (near top of file, imports section):
```groovy
include { SPRING_DECOMPRESS } from './modules/demux'
```

- [ ] **Step 2: Replace the `--from_fastq` block**

Find the current block (main.nf ~lines 341–354):
```groovy
if (params.from_fastq) {
    ch_meta
        .map { meta ->
            def fqs = file("${params.fastq_dir}/**/${meta.id}*.fastq.gz")
            fqs = fqs instanceof List ? fqs : (fqs.exists() ? [fqs] : [])
            def parents = fqs.collect { it.parent.toString() }.unique()
            parents.collect { pdir ->
                def matched_fqs = fqs.findAll { it.parent.toString() == pdir }
                [meta, pdir, matched_fqs]
            }
        }
        .flatMap()
        .filter { meta, fastq_dir, fqs -> !fqs.isEmpty() }
        .set { ch_fastqs }
}
```

Replace with:
```groovy
if (params.from_fastq) {
    // Regular .fastq.gz files — pass through directly
    ch_meta
        .map { meta ->
            def fqs = file("${params.fastq_dir}/**/${meta.id}*.fastq.gz")
            fqs = fqs instanceof List ? fqs : (fqs.exists() ? [fqs] : [])
            def parents = fqs.collect { it.parent.toString() }.unique()
            parents.collect { pdir ->
                def matched = fqs.findAll { it.parent.toString() == pdir }
                [meta, pdir, matched]
            }
        }
        .flatMap()
        .filter { meta, fq_dir, fqs -> !fqs.isEmpty() }
        .set { ch_direct_fastqs }

    // .spring files — decompress first
    ch_meta
        .map { meta ->
            def sprgs = file("${params.fastq_dir}/**/${meta.id}*.spring")
            sprgs = sprgs instanceof List ? sprgs : (sprgs.exists() ? [sprgs] : [])
            sprgs.collect { sp -> [meta, sp] }
        }
        .flatMap()
        .filter { meta, sp -> sp.exists() }
        .set { ch_spring_files }

    SPRING_DECOMPRESS(ch_spring_files)

    // Reformat decompressed output to match ch_fastqs structure: [meta, fq_dir_str, [fqs]]
    SPRING_DECOMPRESS.out.fastqs
        .map { meta, r1, r2 ->
            def fq_dir = r1.parent.toString()
            [meta, fq_dir, [r1, r2]]
        }
        .mix(ch_direct_fastqs)
        .set { ch_fastqs }
}
```

- [ ] **Step 3: Verify compilation**

```bash
nextflow run OSCAR/main.nf --help
# Expected: no groovy errors
```

- [ ] **Step 4: Stub test with a directory containing .spring files**

```bash
nextflow run OSCAR/main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from_fastq true \
    --fastq_dir /path/to/dir/with/spring/files \
    --outdir /tmp/spring_test \
    -stub
# Expected: SPRING_DECOMPRESS appears in DAG; no errors
```

- [ ] **Step 5: Commit**

```bash
git add OSCAR/main.nf
git commit -m "feat: auto-detect and decompress .spring files in --from_fastq mode

When params.fastq_dir contains .spring files matching meta.id, SPRING_DECOMPRESS
runs first to produce .fastq.gz pairs. Mixed with any plain .fastq.gz files found
in the same directory. Both paths produce ch_fastqs in identical format."
```

---

## Task 4: Add SPRING_COMPRESS process to modules/demux.nf

**Files:**
- Modify: `OSCAR/modules/demux.nf` — add process after SPRING_DECOMPRESS

The process:
- Input: `[run_name, r1.fastq.gz, r2.fastq.gz]`
- Output: one `.spring` archive, published to `{outdir}/{run_name}_fastq/`
- Naming: `{r1_prefix}.spring` where prefix = r1 filename with `_R1_001.fastq.gz` stripped

- [ ] **Step 1: Add SPRING_COMPRESS process**

```groovy
// ─── SPRING_COMPRESS ─────────────────────────────────────────────────────────
// Compress a validated R1+R2 FASTQ pair to a single SPRING archive for archival.
// Published to {outdir}/{run_name}_fastq/ replacing the plain FASTQ copy.
// I1/I2 index reads are not compressed (handled by VALIDATE_FASTQ publishDir).

process SPRING_COMPRESS {
    tag "${run_name}/${r1.name}"
    label 'process_medium'
    container "${params.container_spring}"
    publishDir { "${params.outdir}/${run_name}_fastq" }, mode: 'copy'

    input:
    tuple val(run_name), path(r1), path(r2)

    output:
    tuple val(run_name), path("*.spring"), emit: archive
    path "versions.yml",                  emit: versions

    script:
    def prefix = r1.name.replaceAll(/_R1_001\.fastq\.gz$/, '')
    """
    spring -c \\
        -i ${r1} ${r2} \\
        -o ${prefix}.spring \\
        -g \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spring: \$(spring --version 2>&1 | head -1 | sed 's/Spring v//')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add OSCAR/modules/demux.nf
git commit -m "feat: add SPRING_COMPRESS process for archival of BCL→FASTQ output"
```

---

## Task 5: Wire SPRING_COMPRESS into FASTQ_QC subworkflow

**Files:**
- Modify: `OSCAR/subworkflows/fastq_qc.nf`
- Modify: `OSCAR/modules/demux.nf` — change VALIDATE_FASTQ publishDir to skip R1/R2

This task has two parts:
1. VALIDATE_FASTQ stops publishing R1/R2 files (they'll be archived as spring instead)
2. FASTQ_QC subworkflow pairs R1+R2 and calls SPRING_COMPRESS

### Part A: Update VALIDATE_FASTQ publishDir

- [ ] **Step 1: Modify VALIDATE_FASTQ to suppress R1/R2 publishing**

In modules/demux.nf, find the VALIDATE_FASTQ publishDir directive (line ~320):
```groovy
publishDir { "${params.outdir}/${meta.run_name}_fastq" }, mode: 'copy',
    saveAs: { fn -> file(fn).name }
```

Replace with:
```groovy
publishDir { "${params.outdir}/${meta.run_name}_fastq" }, mode: 'copy',
    saveAs: { fn ->
        // R1/R2 reads are archived as .spring by SPRING_COMPRESS; skip plain copy
        file(fn).name =~ /_R[12]_/ ? null : file(fn).name
    }
```

This publishes I1/I2 index reads as-is and suppresses R1/R2 (which get spring-compressed).

### Part B: Add R1/R2 pairing and SPRING_COMPRESS to FASTQ_QC

- [ ] **Step 2: Update subworkflows/fastq_qc.nf**

Add the include and spring compression block. Replace the full file content:

```groovy
include { VALIDATE_FASTQ } from '../modules/demux'
include { FALCO          } from '../modules/demux'
include { SPRING_COMPRESS } from '../modules/demux'

workflow FASTQ_QC {
    take:
        ch_fastqs   // [meta, fastq_dir_str, [fq_files]]

    main:
        // Validate each FASTQ individually as a separate task
        ch_fastqs
            .transpose(by: 2)
            .map { meta, fq_dir, fastq -> [meta, fq_dir, fastq, fastq.name] }
            .set { ch_to_validate }

        VALIDATE_FASTQ(ch_to_validate)

        // Reassemble validated files into per-library lists
        VALIDATE_FASTQ.out.fastq
            .map { meta, fq_dir, fastq -> [meta.id, meta, fq_dir, fastq] }
            .groupTuple(by: 0)
            .map { id, metas, fq_dirs, fastqs ->
                [metas[0], fq_dirs.unique(false), fastqs]
            }
            .set { ch_validated_fastqs }

        // Pair R1+R2 reads for spring compression (I1/I2 excluded by regex)
        VALIDATE_FASTQ.out.fastq
            .filter { meta, fq_dir, fastq -> fastq.name =~ /_R[12]_/ }
            .map { meta, fq_dir, fastq ->
                def pair_key = [meta.run_name, fastq.name.replaceAll(/_R[12]_001\.fastq\.gz$/, '')]
                [pair_key, fastq]
            }
            .groupTuple()
            .map { pair_key, fqs ->
                def fl = []
                fqs.each { fl << it }
                def run_name = pair_key[0]
                def r1 = fl.find { it.name =~ /_R1_/ }
                def r2 = fl.find { it.name =~ /_R2_/ }
                [run_name, r1, r2]
            }
            .filter { run_name, r1, r2 -> r1 && r2 }
            .set { ch_pairs_to_compress }

        SPRING_COMPRESS(ch_pairs_to_compress)

        // FALCO QC on R-reads only (R1/R2/R3; I1/I2 skipped)
        ch_fastqs
            .flatMap { meta, fq_dir, fq_files ->
                def files = fq_files instanceof List ? fq_files : [fq_files]
                files
                    .findAll { f -> f.name =~ /_R[0-9]+_/ && f.size() > 1024 * 1024 }
                    .collect { f -> [meta.run_name, f.name.replaceAll(/\.fastq\.gz$/, ''), f] }
            }
            .set { ch_falco_input }

        FALCO(ch_falco_input)

        FALCO.out.report
            .groupTuple(by: 0)
            .set { ch_falco_reports }

    emit:
        fastqs        = ch_validated_fastqs   // [meta, [fastq_dir_strings], [fastq_files]]
        falco_reports = ch_falco_reports      // [run_name, [report_dirs]]
}
```

- [ ] **Step 3: Verify compilation**

```bash
nextflow run OSCAR/main.nf --help
# Expected: no groovy errors
```

- [ ] **Step 4: Stub run to confirm DAG**

```bash
nextflow run OSCAR/main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir /path/to/bcl \
    --outdir /tmp/spring_compress_test \
    --run_until FASTQ \
    -stub
# Expected: SPRING_COMPRESS appears in DAG; VALIDATE_FASTQ present; no errors
```

- [ ] **Step 5: Commit**

```bash
git add OSCAR/modules/demux.nf OSCAR/subworkflows/fastq_qc.nf
git commit -m "feat: compress BCL→FASTQ R1/R2 pairs to SPRING for archival

VALIDATE_FASTQ publishDir now skips R1/R2 files (saveAs returns null).
FASTQ_QC pairs R1+R2 by sample prefix and calls SPRING_COMPRESS.
I1/I2 index reads still published as .fastq.gz (small, no spring needed).
Downstream (cellranger etc.) unchanged — uses work-dir staged .fastq.gz."
```

---

## Task 6: End-to-end validation

- [ ] **Step 1: Full BCL run with spring compression enabled**

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir /path/to/bcl \
    --outdir /sc-scratch/.../spring_test \
    --run_until FASTQ \
    -resume
# Monitor: check SPRING_COMPRESS tasks complete
```

- [ ] **Step 2: Verify .spring files published**

```bash
ls /sc-scratch/.../spring_test/R???_fastq/*.spring | head -5
# Expected: spring archive per R1+R2 pair
# Example: CITE_CD34_NEU_exp2_libA_S1_L001.spring

# Verify .fastq.gz are NOT published for R1/R2:
ls /sc-scratch/.../spring_test/R???_fastq/*_R1_*.fastq.gz 2>/dev/null
# Expected: no output (empty)

# Verify I-reads still published:
ls /sc-scratch/.../spring_test/R???_fastq/*_I1_*.fastq.gz 2>/dev/null
# Expected: I1 files present (if BCL produced them)
```

- [ ] **Step 3: Verify spring files are valid (spot-check)**

```bash
# Decompress a test spring file and check read count matches original:
apptainer exec /sc-scratch/.../apptainer_cache/quay.io-biocontainers-spring_1.1.1--h4ac6f70_3.img \
    spring -d -i test_S1_L001.spring \
           -o /tmp/check_R1.fastq.gz /tmp/check_R2.fastq.gz -g
zcat /tmp/check_R1.fastq.gz | wc -l
# Lines / 4 = read count; should match the validated FASTQ
```

- [ ] **Step 4: Validate --from_fastq with spring input**

```bash
# Replace the published fastq dir with a spring-only dir and re-run:
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from_fastq true \
    --fastq_dir /sc-scratch/.../spring_test/R???_fastq \
    --outdir /sc-scratch/.../spring_round_trip_test \
    --run_until cellranger \
    -resume
# Expected: SPRING_DECOMPRESS runs for each .spring file; cellranger receives .fastq.gz
```

---

## File Change Summary

| File | Change |
|---|---|
| `nextflow.config` | Add `container_spring` |
| `modules/demux.nf` | Add `SPRING_DECOMPRESS` process; add `SPRING_COMPRESS` process; modify `VALIDATE_FASTQ` publishDir `saveAs` to skip R1/R2 |
| `subworkflows/fastq_qc.nf` | Add `SPRING_COMPRESS` include + R1/R2 pairing channel + `SPRING_COMPRESS` call |
| `main.nf` | Add `SPRING_DECOMPRESS` include; replace `--from_fastq` block with spring-aware version |

## Key design invariant

Downstream of `ch_fastqs`, **nothing changes**. FASTQ_QC still emits `[meta, [fq_dirs], [fqs]]` with plain `.fastq.gz` files staged in work dirs. Cellranger, CellBender, and all other processes receive the same channel structure as before. Spring only affects (a) what gets published for archival and (b) what gets loaded when `--from_fastq` is used.
