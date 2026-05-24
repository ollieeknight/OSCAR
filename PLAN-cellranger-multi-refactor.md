# Plan: CELLRANGER_MULTI Refactor

Fixes three bugs in the OSCAR pipeline's GEX counting path:
1. Missing `groupTuple` causing wrong types in channel (root cause of ArrayBag errors)
2. ArrayBag-unsafe materialisation in the process script block
3. Fragile heredoc CSV generation → Python template

---

## Phase 0: Discovery Summary (Complete)

### Key files

| File | Lines | Role |
|---|---|---|
| `main.nf` | 346–387 | Channel construction → COUNT_GEX |
| `modules/count_gex.nf` | 1–87 | CELLRANGER_MULTI process |
| `subworkflows/count_gex.nf` | 1–13 | Thin wrapper, no logic |
| `subworkflows/demux.nf` | 33–43 | Emits `[meta, fq_dir_str, [fqs]]` per sample |

### Root cause: missing groupTuple (main.nf:359)

`ch_routed.gex` inherits the branch shape from `ch_fastqs`, which is
`[meta, fastq_dir_string, [fqs]]` — one tuple per sample row.

The existing `.map { lid, metas, fastq_dirs ->` at line 360 assumes the channel
is already grouped (`lid` = library_id string, `metas` = list, `fastq_dirs` = list).
Without the grouping step this binds:
- `lid`       → a single `meta` Map  (wrong type)
- `metas`     → a single fastq_dir String  (wrong type)
- `fastq_dirs`→ a single `[fqs]` List  (wrong type)

The ATAC branch immediately below (lines 380–383) shows the correct idiom:

```groovy
ch_routed.atac
    .map { meta, fastq_dir, fqs -> [meta.library_id, meta, fastq_dir] }
    .groupTuple(by: 0)
    .map { lid, metas, fastq_dirs -> [metas[0], fastq_dirs.toSet().toList()] }
```

The GEX path is missing the first two lines of that pattern.

### Secondary cause: ArrayBag in process (count_gex.nf:16, 49)

After the groupTuple is restored, Nextflow passes `metas` and `fastq_dirs` as
`nextflow.util.ArrayBag`. These lines then throw `UnsupportedOperationException`:

```groovy
// line 16 — [metas].flatten() on an ArrayBag fails
def metas_list = [metas].flatten().findAll { it != null }

// line 49 — `as ArrayList` cast on an ArrayBag fails
def dirs_list = (fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]) as ArrayList
```

### Tertiary cause: heredoc CSV fragility (count_gex.nf:66–85)

GString interpolation of multi-line section strings (`ge_section`, `vdj_section`,
`feature_section`) into a heredoc produces unpredictable blank lines when optional
sections are empty strings. The CSV parser in `cellranger multi` will error on blank
lines between sections.

### Allowed Nextflow operators (from CLAUDE.md)

`.branch`, `.join`, `.groupTuple`, `.combine`, `.flatten`, `.multiMap`, `.mix`,
`.collect`, `.map`, `.filter`, `.set`, `groupKey()`

### Safe materialisation pattern (from main.nf:363–366, already correct)

```groovy
def unique_metas = []
metas.each { m -> if (seen.add(m.modality)) unique_metas << m }
```

Use `<<` append into a literal `[]` — never `findAll`, `as ArrayList`, or
`.flatten()` on a value that may be an ArrayBag.

---

## Phase 1: Fix main.nf — Add Missing groupTuple

**File**: `main.nf`  
**Scope**: lines 359–360 only (insert two lines before the existing map)

### What to change

Current (line 359):
```groovy
            ch_routed.gex
                .map { lid, metas, fastq_dirs ->
```

Replace with:
```groovy
            ch_routed.gex
                .map { meta, fastq_dir, fqs -> [meta.library_id, meta, fastq_dir] }
                .groupTuple(by: 0)
                .map { lid, metas, fastq_dirs ->
```

The body of the inner `.map` (lines 361–376) already uses safe `[]`-append
patterns for `unique_metas` and `unique_dirs` and will work correctly once
`metas` and `fastq_dirs` are proper ArrayBag inputs.

### What NOT to change

- Do not touch lines 361–376 (dedup body) — already correct
- Do not change the ATAC or ASAP branches
- Do not add `groupKey()` — not needed here (no size hint required)

### Verification checklist

```bash
# Confirm two groupTuple calls for gex (new) and atac (existing)
grep -n 'groupTuple\|ch_routed\.gex' main.nf

# Stub dry-run — validates channel graph compiles (no FASTQ needed)
nextflow run main.nf -profile standard --from_fastq true \
    --samplesheet assets/example_metadata.csv \
    --fastq_dir /dev/null --outdir /tmp/oscar_stub_test -stub 2>&1 | tail -30
```

Expected: `groupTuple` appears at two nearby lines for `ch_routed.gex`, no
`WARN` or `ERROR` in stub output.

---

## Phase 2: Fix count_gex.nf — ArrayBag-Safe Materialisation

**File**: `modules/count_gex.nf`  
**Scope**: lines 16 and 49

### What to change

**Line 16** — replace:
```groovy
    def metas_list   = [metas].flatten().findAll { it != null }
    if (!metas_list) error "CELLRANGER_MULTI: metas is empty for library ${library_id}"
```
with:
```groovy
    def metas_list = []
    metas.each { m -> if (m != null) metas_list << m }
    if (!metas_list) error "CELLRANGER_MULTI: metas is empty for library ${library_id}"
```

**Line 49** — replace:
```groovy
    def dirs_list = (fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]) as ArrayList
```
with:
```groovy
    def dirs_list = []
    fastq_dirs.each { d -> dirs_list << d }
    if (!dirs_list) error "CELLRANGER_MULTI: fastq_dirs is empty for library ${library_id}"
```

All subsequent collection ops (lines 25–27, 51–63) operate on `metas_list`
and `dirs_list` — both plain `java.util.ArrayList` after the fix — so no
further changes are needed in the script block.

### What NOT to change

- Do not touch the section-building logic (lines 31–64)
- Do not change the process input/output signatures
- Do not add defensive copies inside the section-builder vars

### Verification checklist

```bash
# Zero remaining ArrayBag-unsafe patterns
grep -n 'flatten\|as ArrayList\|ArrayBag' modules/count_gex.nf
# Expected: 0 results

# Confirm metas_list and dirs_list are the only collection vars used downstream
grep -n 'metas_list\|dirs_list\|metas\b\|fastq_dirs\b' modules/count_gex.nf
# Expected: metas and fastq_dirs only appear in input: and in the two
# materialisation blocks; all other uses reference metas_list / dirs_list
```

---

## Phase 3: Python Config Template

**New file**: `bin/cellranger_multi_config.py`  
**Modified file**: `modules/count_gex.nf` (script: block, lines 66–85)

### Context

The nf-core scrnaseq pipeline uses `template "cellranger_multi.py"` to write
the CSV. OSCAR cannot use Nextflow templates (no `templates/` in module path)
but can use an executable script in `bin/` — Nextflow automatically adds `bin/`
to `PATH` inside every process work dir.

### bin/cellranger_multi_config.py

Write a Python 3 script (stdlib only — must work inside the cellranger container)
that accepts CLI args and writes `multi_config.csv`.

**CLI signature**:
```
cellranger_multi_config.py \
    --ref_gex    PATH \
    --create_bam true|false \
    --modalities GEX[,ADT][,HTO][,VDJ-T][,VDJ-B][,CRISPR] \
    --assay      CITE|DOGMA|Flex|Multiome|GEX|ASAP \
    --chemistry  SC3Pv3|ARC-v1|... \
    --adt_csv    PATH|NO_FILE \
    --ref_vdj    PATH \
    --output     multi_config.csv
```

**Section logic** (mirrors current Groovy exactly):

```
[gene-expression]
reference,{ref_gex}
create-bam,{create_bam}
chemistry,ARC-v1          ← only if assay in {DOGMA, Multiome}
chemistry,{chemistry}     ← only if assay == Flex
                          ← omit chemistry line otherwise

[vdj]                     ← only if VDJ-T or VDJ-B in modalities
reference,{ref_vdj}

[feature]                 ← only if (ADT or HTO in modalities) AND adt_csv != NO_FILE
reference,{adt_csv}

[libraries]
fastq_id,fastqs,feature_types
                          ← rows appended by bash lib_check at runtime
```

The script must:
- Exit 0 on success (CSV written to `--output` path)
- Exit 1 with a descriptive message if required args are missing
- Write NO blank lines between sections (strip trailing newlines from each section)
- Be executable (`#!/usr/bin/env python3`, `chmod +x`)

**Minimal implementation skeleton**:
```python
#!/usr/bin/env python3
import argparse, sys

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--ref_gex',    required=True)
    p.add_argument('--create_bam', required=True)
    p.add_argument('--modalities', required=True)
    p.add_argument('--assay',      required=True)
    p.add_argument('--chemistry',  default='')
    p.add_argument('--adt_csv',    required=True)
    p.add_argument('--ref_vdj',    required=True)
    p.add_argument('--output',     required=True)
    args = p.parse_args()

    mods = [m.strip() for m in args.modalities.split(',')]
    is_dogma_or_multiome = args.assay in ('DOGMA', 'Multiome')
    is_flex              = args.assay == 'Flex'
    has_vdj    = any(m in mods for m in ('VDJ-T', 'VDJ-B'))
    has_feature = any(m in mods for m in ('ADT', 'HTO')) and args.adt_csv != 'NO_FILE'

    lines = ['[gene-expression]',
             f'reference,{args.ref_gex}',
             f'create-bam,{args.create_bam}']
    if is_dogma_or_multiome:
        lines.append('chemistry,ARC-v1')
    elif is_flex:
        lines.append(f'chemistry,{args.chemistry}')

    if has_vdj:
        lines += ['', '[vdj]', f'reference,{args.ref_vdj}']

    if has_feature:
        lines += ['', '[feature]', f'reference,{args.adt_csv}']

    lines += ['', '[libraries]', 'fastq_id,fastqs,feature_types']

    with open(args.output, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
```

Note: blank lines between sections ARE correct in the cellranger multi CSV format
(section delimiter). The bug in the heredoc was about *extra* or *missing* blank
lines caused by GString empty-section interpolation. The Python version controls
blank-line placement explicitly.

### modules/count_gex.nf — script: block update

Replace lines 66–85 (heredoc + cellranger invocation) with:

```bash
    cellranger_multi_config.py \
        --ref_gex    "${ref_gex}" \
        --create_bam "${create_bam}" \
        --modalities "${metas_list.collect { it.modality }.join(',')}" \
        --assay      "${meta.assay}" \
        --chemistry  "${meta.chemistry ?: ''}" \
        --adt_csv    "${adt_csv.toAbsolutePath()}" \
        --ref_vdj    "${ref_vdj}" \
        --output     multi_config.csv

    ${lib_check_script}

    cellranger multi \\
        --id        "${library_id}" \\
        --csv       multi_config.csv \\
        --localcores ${task.cpus} \\
        --localmem  ${task.memory.toGiga()}

    rm -rf "${library_id}/SC_MULTI_CS" "${library_id}/_"*

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | head -1 | sed 's/.* //')
END_VERSIONS
```

Keep `lib_check_script` (lines 50–64, Groovy) — it generates bash that evaluates
actual file sizes at runtime, which cannot be moved to Python.

### What NOT to change

- Do not change the process `input:` / `output:` signatures
- Do not move `lib_check_script` generation out of Groovy — it needs runtime bash eval
- Do not use Nextflow `template` directive — OSCAR has no `templates/` directory
- Do not add new params for ref paths — already available as `params.*` inside process

### Verification checklist

```bash
# Script is executable and parses args
python3 bin/cellranger_multi_config.py --help

# GEX-only (no VDJ, no ADT)
python3 bin/cellranger_multi_config.py \
    --ref_gex /ref/human --create_bam false \
    --modalities GEX --assay GEX --chemistry SC3Pv3 \
    --adt_csv NO_FILE --ref_vdj /ref/vdj \
    --output /tmp/test_gex.csv && cat /tmp/test_gex.csv

# CITE (GEX + ADT)
python3 bin/cellranger_multi_config.py \
    --ref_gex /ref/human --create_bam false \
    --modalities GEX,ADT --assay CITE --chemistry SC3Pv3 \
    --adt_csv /path/to/adt.csv --ref_vdj /ref/vdj \
    --output /tmp/test_cite.csv && cat /tmp/test_cite.csv

# DOGMA (GEX + ADT, chemistry=ARC-v1)
python3 bin/cellranger_multi_config.py \
    --ref_gex /ref/human --create_bam false \
    --modalities GEX,ADT --assay DOGMA --chemistry ARC-v1 \
    --adt_csv /path/to/adt.csv --ref_vdj /ref/vdj \
    --output /tmp/test_dogma.csv && cat /tmp/test_dogma.csv

# Flex
python3 bin/cellranger_multi_config.py \
    --ref_gex /ref/human --create_bam false \
    --modalities GEX --assay Flex --chemistry SC-FB-v2 \
    --adt_csv NO_FILE --ref_vdj /ref/vdj \
    --output /tmp/test_flex.csv && cat /tmp/test_flex.csv

# Confirm no heredoc in updated process
grep -n 'MULTIEOF\|<< .MULTI\|cat > multi_config' modules/count_gex.nf
# Expected: 0 results
```

---

## Phase 4 (Optional): Improve lib_check_script

**File**: `modules/count_gex.nf` lines 58–61  
**Decision**: Apply if GNU `find -printf` portability is a concern; skip otherwise.  
OSCAR runs in Linux containers so the current approach works — this is cosmetic.

Current (GNU-only, computes size sum):
```bash
if [ $(find "${dir}" -maxdepth 2 -name "${m.id}*.fastq.gz" \
    -printf '%s\n' 2>/dev/null | awk '{s+=$1} END{printf "%.0f\n", s}') -ge 10485760 ]; then
    echo "${m.id},${dir},${ft}" >> multi_config.csv
fi
```

Simpler POSIX equivalent (existence check, no awk):
```bash
if find "${dir}" -maxdepth 2 -name "${m.id}*.fastq.gz" \
    -print -quit 2>/dev/null | grep -q .; then
    echo "${m.id},${dir},${ft}" >> multi_config.csv
fi
```

Same semantics (modality present ↔ at least one matching FASTQ exists), faster,
no GNU dependency, no size threshold magic number.

---

## Phase 5: Final Verification

Run after all phases are complete.

```bash
# 1. Channel graph
grep -n 'groupTuple' main.nf
# Expected: 3 lines — gex (new), atac, asap_adt

# 2. No ArrayBag-unsafe patterns in process
grep -n 'flatten\|as ArrayList\|ArrayBag' modules/count_gex.nf
# Expected: 0 results

# 3. No heredoc CSV in process
grep -n 'MULTIEOF\|cat > multi_config' modules/count_gex.nf
# Expected: 0 results

# 4. Python script exists and is executable
ls -la bin/cellranger_multi_config.py

# 5. Stub dry-run compiles cleanly
nextflow run main.nf -profile standard \
    --from_fastq true \
    --samplesheet assets/example_metadata.csv \
    --fastq_dir /dev/null \
    --outdir /tmp/oscar_stub \
    -stub 2>&1 | grep -E 'ERROR|WARN|error'
# Expected: 0 lines

# 6. On first real run: inspect generated CSV
cat work/<hash>/multi_config.csv
# Expected: correct sections, one blank line between sections,
#           no double blank lines, no trailing blank lines before [libraries]
```

---

## Execution Order

Phases are independent within each phase boundary but must be executed in order:

```
Phase 1 (main.nf groupTuple)
    → Phase 2 (count_gex.nf ArrayBag fix)       [depends on Phase 1 to understand what metas is]
    → Phase 3 (Python template)                  [depends on Phase 2 — process already fixed]
    → Phase 4 (optional lib_check improvement)   [independent of 1–3]
    → Phase 5 (verification)                     [requires all prior phases]
```

Phase 1 alone fixes the wrong-type channel bug. Phase 2 alone does nothing
(no groupTuple = no ArrayBag to defend against). Both must be applied together
for the pipeline to run correctly.
