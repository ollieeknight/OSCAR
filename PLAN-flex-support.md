# PLAN: Full Flex (Fixed RNA Profiling) Support

**Goal**: Complete Flex support in OSCAR — probe set resolution, multiplexed sample barcode
assignments, correct config generation, validation, and custom probe set option.

**Constraint**: Cannot add samplesheet columns (must match `bash/metadata.csv`).

---

## Phase 0: Documentation Discovery (DONE)

### 10x Genomics cellranger multi: Flex requirements

Source: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts

**`[gene-expression]` section — Flex-specific:**
```
probe-set,/abs/path/to/probe-set.csv     ← REQUIRED for all Flex
filter-probes,true                        ← optional; default true (exclude deprecated probes)
reference,/abs/path/to/transcriptome     ← optional for Flex v2, required for v1
create-bam,false                          ← recommended false for Flex
chemistry,Flex-v2-R1                      ← already implemented
```

Built-in probe set files (inside `cellranger-x.y.z/probe_sets/`):
- Flex v2 human: `Chromium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv`
- Flex v2 mouse: `Chromium_Mouse_Transcriptome_Probe_Set_v2.0_GRCm39-2024-A.csv`
- Flex v1 human: `Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv`
- Flex v1 mouse: `Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv`

**`[samples]` section — multiplexed Flex only:**
```
[samples]
sample_id,probe_barcode_ids
Sample1,A-A01              ← Flex v2 format: <plate>-<well>
Sample2,A-B01
Sample3,A-A01|A-A02        ← pipe-separated for multiple barcodes per sample
```
- Flex v1 format: `BCxxx` (e.g. `BC001`); with ADT: `BC001+AB005`
- Flex v2 format: `<plate>-<well>` (e.g. `A-A01`)
- Singleplex (n_donors=1): `[samples]` section is NOT required

### Current pipeline state (OSCAR)

**Working:**
- Flex assay routes to GEX pipeline (main.nf:336, 429)
- Demux OverrideCycles hardcoded to `DI_Flex-v2_GEX` (modules/demux.nf:73)
- Chemistry line added: `chemistry,Flex-v2-R1` (main.nf:482-483)
- RNA velocity skipped for Flex (probe-based, no intronic signal)

**Missing (critical):**
1. No `probe-set` line in `[gene-expression]` section → cellranger multi will fail
2. No `[samples]` section for multiplexed Flex (n_donors > 1)
3. No probe set params in nextflow.config
4. `create-bam` not overridden to `false` for Flex
5. Flex-specific validation (index_type must be DI, modality must be GEX)

### Design decision: probe barcode assignments without new samplesheet columns

Cannot add columns. Solution: **reuse `adt_file`** column when `assay=Flex`.
- When `assay=Flex` and `adt_file` is set → resolve as a **flex samples CSV** (not a feature barcode CSV)
- Flex samples CSV format: `sample_id,probe_barcode_ids[,description]`
- Same 3-tier resolution logic as ADT CSV (local → parent → centralized fallback)
- When `assay=Flex` and `adt_file` is empty → singleplex mode (no `[samples]` section)
- Existing Flex samplesheet entries (`adt_file` blank) remain valid with no changes

This reuse is safe because for Flex, `adt_file` has always been blank (Flex has no antibody
capture). Any existing Flex entries will continue to work unchanged.

---

## Phase 1: nextflow.config — Probe Set Parameters

**File**: `nextflow.config`

Add probe set reference paths to `params {}` block, after the existing `ref_*` params (line 29):

```groovy
// Flex probe sets (inside the cellranger installation)
probeset_human_v2      = "/home/knighto/bin/cellranger-10.0.0/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv"
probeset_mouse_v2      = "/home/knighto/bin/cellranger-10.0.0/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v2.0_GRCm39-2024-A.csv"
probeset_human_v1      = "/home/knighto/bin/cellranger-10.0.0/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"
probeset_mouse_v1      = "/home/knighto/bin/cellranger-10.0.0/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv"
probeset_custom        = null   // override: absolute path to custom probe set CSV
```

**Probe set auto-selection logic** (to implement in main.nf helper function):
- `probeset_custom` set → use that (all species, all chemistry variants)
- else chemistry starts with `Flex-v2` → v2 probe set keyed by species
- else (Flex-v1 / SFRP / MFRP) → v1 probe set keyed by species

**Verification**: `grep -n 'probeset' nextflow.config` returns 5 lines.

---

## Phase 2: main.nf — Helper Function + Meta Map Extension

**File**: `main.nf`

### 2a. Add `get_flex_probeset(meta)` helper function

After `get_viral_whitelist()` (around line 183), add:

```groovy
def get_flex_probeset(meta) {
    if (params.probeset_custom)
        return params.probeset_custom
    def is_v2 = meta.chemistry?.startsWith('Flex-v2')
    if (meta.species == 'human')
        return is_v2 ? params.probeset_human_v2 : params.probeset_human_v1
    if (meta.species == 'mouse')
        return is_v2 ? params.probeset_mouse_v2 : params.probeset_mouse_v1
    error "get_flex_probeset: unknown species '${meta.species}'"
}
```

### 2b. Extend `parse_row()` to handle Flex samples file

In `parse_row()` (main.nf:122-166), the existing adt_csv resolution block (lines 131-147)
already resolves `adt_file` → `adt_csv_path`. That logic stays unchanged.

Add to the returned map (after `adt_csv_path`):

```groovy
flex_samples_path: (row.assay.trim().equalsIgnoreCase('Flex') && adt_csv_path) ? adt_csv_path : null,
adt_csv_path:      (row.assay.trim().equalsIgnoreCase('Flex')) ? null : adt_csv_path
```

**Rationale**: For Flex, `adt_csv_path` is redirected to `flex_samples_path`. The
`[feature]` section check (`has_adt && has_adt_csv`) already gates on ADT/HTO modality
being present, so setting `adt_csv_path=null` for Flex is safe — Flex has no ADT modality.

Also, update the existing adt_csv warning message to conditionally mention Flex:
```groovy
def warn_context = row.assay.trim().equalsIgnoreCase('Flex') ? 'flex_samples' : 'ADT'
log.warn "WARNING: ${warn_context} file '${adt_file}.csv' not found ..."
```

**Verification**: 
- `grep -n 'flex_samples_path' main.nf` — appears in parse_row return map
- `grep -n 'adt_csv_path' main.nf` — still present in non-Flex path

---

## Phase 3: main.nf — Config Header Generation

**File**: `main.nf` (config building block, lines ~479-494)

This is where the `[gene-expression]`, `[vdj]`, `[feature]` sections are built as a
`List<String>`. All changes happen here in the `.map {}` block.

### 3a. Add `probe-set` line for Flex

After the `chemistry` line addition (line 483), add:

```groovy
if (meta.assay == 'Flex') {
    def ps = get_flex_probeset(meta)
    if (!ps) error "Flex library '${lid}': probe-set not resolved. " +
                   "Set params.probeset_human_v2/probeset_mouse_v2 or params.probeset_custom."
    lines << "probe-set,${ps}"
}
```

### 3b. Override `create-bam` to `false` for Flex

The config header currently adds `create-bam,${create_bam}`. Wrap in conditional:

```groovy
def bam_flag = (meta.assay == 'Flex') ? 'false' : create_bam
def lines = ['[gene-expression]',
             "reference,${ref_gex}",
             "create-bam,${bam_flag}"]
```

**Rationale**: 10x docs explicitly recommend `create-bam,false` for Flex. BAM generation
for Flex can exceed 1TB and is rarely needed for typical GEX analysis.

### 3c. Add `[samples]` section from flex_samples_path

After the `[feature]` section block, add:

```groovy
// Flex multiplexed: add [samples] section from flex samples file
def flex_samples_path = all_metas.collect { it.flex_samples_path }.find { it }
if (meta.assay == 'Flex' && flex_samples_path) {
    def samples_csv = new File(flex_samples_path)
    if (!samples_csv.exists())
        error "Flex library '${lid}': flex_samples file not found: ${flex_samples_path}"
    def samples_headers = samples_csv.readLines()[0]
    lines << ''
    lines << '[samples]'
    lines << samples_headers
    samples_csv.readLines().tail().each { row_line ->
        if (row_line.trim()) lines << row_line
    }
} else if (meta.assay == 'Flex' && meta.n_donors > 1) {
    log.warn "WARNING: Flex library '${lid}' has n_donors=${meta.n_donors} but no flex_samples " +
             "file provided (adt_file column is empty). Cellranger will treat as singleplex. " +
             "Provide a flex_samples CSV via the adt_file column to enable demultiplexing."
}
```

**Important**: Follow the existing ArrayBag materialisation pattern before iterating
`all_metas`. Use:
```groovy
def ml = []
all_metas.each { m -> ml << m }
def flex_samples_path = ml.collect { it.flex_samples_path }.find { it }
```

**Verification**:
- Run: `grep -A5 '\[samples\]' results/*/outs/multi_config.csv` on a multiplexed Flex run
- Config should have: `[gene-expression]`, `probe-set`, `create-bam,false`, `chemistry`,
  `[libraries]`, `[samples]`

---

## Phase 4: Validation — preflight_samplesheet()

**File**: `main.nf` (lines 84-120)

Add Flex-specific validation inside the `lines.tail().eachWithIndex` loop, after the
existing chemistry check (line 119):

```groovy
// Flex-specific constraints
if (row.assay.equalsIgnoreCase('Flex')) {
    if (row.index_type?.trim() != 'DI')
        error "ERROR: row ${i + 2}: Flex assay requires index_type=DI, got '${row.index_type}'"
    if (row.modality?.trim() != 'GEX')
        error "ERROR: row ${i + 2}: Flex assay only supports modality=GEX, got '${row.modality}'"
}
```

**Rationale**: Demux already hardcodes DI for Flex (modules/demux.nf:73) — enforcing this
at validation gives a clear error instead of silent override.

**Verification**: Samplesheet with `Flex,...,SI,...` should fail with clear error message.

---

## Phase 5: Documentation — CLAUDE.md

**File**: `CLAUDE.md`

### 5a. Update samplesheet examples section

Add/update the Flex example:

```csv
# Flex singleplex (1 sample per GEM well, no demux)
Flex,IBD,1,A,GEX,Flex-v2-R1,DI,SI-GA-B1,human,1,

# Flex multiplexed (4 samples per GEM well — adt_file points to flex_samples CSV)
Flex,IBD,2,A,GEX,Flex-v2-R1,DI,SI-GA-B1,human,4,IBD_flex_samples
```

Flex samples CSV format (placed at `{samplesheet_dir}/adt_files/IBD_flex_samples.csv`):
```csv
sample_id,probe_barcode_ids,description
IBD_Ctrl1,A-A01,Control replicate 1
IBD_Ctrl2,A-B01,Control replicate 2
IBD_Treat1,A-C01,Treated replicate 1
IBD_Treat2,A-D01,Treated replicate 2
```

Probe barcode ID formats:
- Flex v2: `<plate>-<well>` (e.g. `A-A01`); pipe-separated for sub-pooled: `A-A01|A-A02`
- Flex v1: `BCxxx` (e.g. `BC001`); with ADT: `BC001+AB005`

### 5b. Update samplesheet column table

Add note to `adt_file` row:
```
adt_file | stem of ADT CSV for CITE/DOGMA; OR stem of flex_samples CSV when assay=Flex
         | (see ADT CSV resolution for lookup logic); blank = singleplex Flex or no ADT
```

### 5c. Update [gene-expression] config section docs

Update the config section table to include:
| Assay | Added lines |
|-------|-------------|
| Flex  | `probe-set,{auto-resolved}`, `chemistry,{Flex-v2-R1|...}`, `create-bam,false` |

### 5d. Update Key params section

Add:
```groovy
// ── Flex probe sets ────────────────────────────────────────────────────────────
probeset_human_v2 = "/home/knighto/bin/cellranger-10.0.0/probe_sets/..."
probeset_mouse_v2 = "/home/knighto/bin/cellranger-10.0.0/probe_sets/..."
probeset_human_v1 = "/home/knighto/bin/cellranger-10.0.0/probe_sets/..."
probeset_mouse_v1 = "/home/knighto/bin/cellranger-10.0.0/probe_sets/..."
probeset_custom   = null   // set to override auto-selection
```

### 5e. Update What NOT to do section

Add:
- Do not add `probe-set` to the global params samplesheet — it is resolved automatically from species + chemistry
- Flex `adt_file` column is reused for flex_samples CSV (probe barcode assignments), NOT for antibody capture

---

## Phase 6: Verification Checklist

### Static checks (grep)
```bash
cd OSCAR
grep -n 'probe-set\|probeset\|flex_samples' main.nf
grep -n 'probeset' nextflow.config
grep -n 'Flex' modules/demux.nf
```

### Config generation dry-run

1. **Singleplex Flex** (no adt_file):
   - Expected config: `[gene-expression]` with `reference`, `create-bam,false`,
     `chemistry,Flex-v2-R1`, `probe-set,/path/to/v2_human_probeset.csv`;
     `[libraries]` only; **no** `[samples]` section.

2. **Multiplexed Flex** (adt_file pointing to flex_samples CSV):
   - Expected: same as above + `[samples]` section with sample_id/probe_barcode_ids.

3. **Custom probe set** (`--probeset_custom /path/to/custom.csv`):
   - Expected: `probe-set,/path/to/custom.csv` regardless of species.

### Validation checks

4. Samplesheet with `Flex,...,SI,...` → error: "Flex assay requires index_type=DI"
5. Samplesheet with `Flex,...,ADT,...` → error: "Flex assay only supports modality=GEX"
6. Samplesheet with `Flex,...,n_donors=4` + empty adt_file → warning logged, no error

### Backward compatibility

7. Existing CITE/DOGMA entries with `adt_file` set → `adt_csv_path` still resolved
   correctly; `flex_samples_path` is null (assay != Flex)
8. Existing Flex entry with empty `adt_file` (from CLAUDE.md example) → singleplex mode,
   probe-set auto-resolved from species+chemistry, run succeeds

---

## Anti-patterns to Avoid

- Do NOT call `.collect {}` or `.any {}` on ArrayBag directly in a process `script:` block
  — iterate into a `def ml = []` in the `.map {}` operator first
- Do NOT pass the flex_samples CSV as a `path()` input to MULTI_CONFIG — Groovy reads it
  in the `.map {}` block before the process, consistent with how adt_csv is handled
- Do NOT make `probe-set` line conditional on some flag — it is REQUIRED for all Flex runs;
  fail loudly if unresolvable
- Do NOT hardcode probe set filenames — always derive from `params.*` so the cellranger
  version can be updated independently
- Do NOT add `filter-probes,false` by default — leave it to cellranger's default (`true`)
  unless the user explicitly needs deprecated probes
