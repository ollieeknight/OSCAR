# PLAN: RNA Velocity (Spliced/Unspliced Counts) — OSCAR Pipeline

## Settled Decisions

| Decision | Choice | Reason |
|---|---|---|
| Tool | simpleaf USA mode | Best accuracy; already in pipeline via VIRAL_DETECT |
| Input | GEX FASTQs (re-align) | More accurate than BAM-based tools |
| Whitelist | Cellbender barcodes | Pre-filtered cells; smaller output; no 3M file needed |
| Run order | After cellbender | Needs cellbender barcodes as permitted list |
| Output location | `{library_id}/velocity/` | Matches cellbender/, vireo/ sibling pattern |
| R import | `fishpond::loadFry()` | Seurat-compatible; handles USA mode natively |
| Excluded assays | Flex | Probe-based; no intronic signal |
| t2g file | Derived from index dir | `t2g_3col.tsv` is inside `index/` — no separate param |

---

## Index Status

**Human index: BUILT**
```
/home/knighto/work/ref/hs/GRCh38-hardmasked-optimised-arc-simpleaf/
├── index/
│   ├── piscem_idx.*       ← piscem index files
│   ├── simpleaf_index.json
│   └── t2g_3col.tsv       ← t2g lives here; no separate param needed
├── ref/
│   └── roers_ref.fa
└── simpleaf_index_log.json
```
Param: `spliceu_index_human = "/home/knighto/work/ref/hs/GRCh38-hardmasked-optimised-arc-simpleaf/index"`

**Mouse index: NOT YET BUILT** — run Phase 0 below.

---

## Phase 0: Build Mouse spliceu Index (Manual, One-Time)

```bash
apptainer exec $APPTAINER_BIND \
    ~/scratch/apptainer_cache/quay.io-biocontainers-simpleaf-0.25.0--hd612981_0.img \
    simpleaf index \
        --output /home/knighto/work/ref/mm/GRCm38-hardmasked-optimised-arc-simpleaf/ \
        --fasta  /sc-projects/sc-proj-cc12-ag-romagnani/ref/mm/GRCm38-hardmasked-optimised-arc/fasta/genome.fa \
        --gtf    /sc-projects/sc-proj-cc12-ag-romagnani/ref/mm/GRCm38-hardmasked-optimised-arc/genes/genes.gtf.gz \
        --ref-type spliced+unspliced \
        --threads $(nproc)
```

Param after: `spliceu_index_mouse = "/home/knighto/work/ref/mm/GRCm38-hardmasked-optimised-arc-simpleaf/index"`

**Verify:**
- [ ] `ls /home/knighto/work/ref/mm/GRCm38-hardmasked-optimised-arc-simpleaf/index/t2g_3col.tsv` exists
- [ ] `ls /home/knighto/work/ref/mm/GRCm38-hardmasked-optimised-arc-simpleaf/index/piscem_idx.ctab` exists

---

## Allowed simpleaf APIs

Verified from VIRAL_DETECT (`modules/qc.nf:265-274`) and the index run above.

```bash
export ALEVIN_FRY_HOME=${PWD}/.alevin_fry_home
mkdir -p "${ALEVIN_FRY_HOME}"
simpleaf set-paths

simpleaf quant \
    --reads1        "r1a.fastq.gz,r1b.fastq.gz" \
    --reads2        "r2a.fastq.gz,r2b.fastq.gz" \
    --threads       ${task.cpus} \
    --index         /path/to/index \
    --chemistry     10xv3 \
    --t2g-map       /path/to/index/t2g_3col.tsv \
    --resolution    cr-like \
    --unfiltered-pl /path/to/barcodes.txt \
    --usa-mode \
    --output        velocity/
```

**`--unfiltered-pl` with cellbender barcodes**: simpleaf only counts barcodes present in this list. Cellbender CSV has barcodes with `-1` suffix (CellRanger format) — must strip before passing:
```bash
sed 's/-1$//' output_cell_barcodes.csv > barcodes_clean.txt
```

**Chemistry mapping** (meta.chemistry → simpleaf --chemistry):
| OSCAR chemistry | simpleaf | Notes |
|---|---|---|
| SC3Pv2 | 10xv2 | |
| SC3Pv3 | 10xv3 | Most common |
| SC3Pv4 | 10xv4-3p | |
| SC5P, SC5Pv3 | 10xv3 | Same barcode structure as 3' v3 |
| ARCv1 | 10xv3 | Multiome GEX uses v3 barcodes |
| Flex-v2-R1 | — | SKIP; probe-based |

---

## Phase 1: `SIMPLEAF_VELOCITY` Process

**File**: `modules/qc.nf`  
**Location**: After VIRAL_DETECT process (line ~282).  
**Pattern**: MGATK2 for publishDir (writes named dir, publishes at library level); VIRAL_DETECT for simpleaf quant.

```nextflow
// ─── SIMPLEAF_VELOCITY ────────────────────────────────────────────────────────
// Spliced/unspliced quantification for RNA velocity.
// Re-quantifies GEX FASTQs using simpleaf USA mode with a spliceu reference.
// Runs after cellbender; uses cellbender cell barcodes as the permitted list.
// Not run for Flex (probe-based chemistry, no intronic signal).
// R import: fishpond::loadFry("velocity/af_quant", outputFormat = "velocity")

process SIMPLEAF_VELOCITY {
    tag "$meta.library_id"
    label 'process_high'
    container "${params.container_simpleaf}"
    publishDir { "${params.outdir}/${meta.run_name}_outs/${meta.library_id}" },
               mode: 'copy'

    input:
    tuple val(meta), val(gex_fastq_dirs), val(simpleaf_chemistry), path(barcodes)
    path spliceu_index

    output:
    tuple val(meta), path("velocity/"), emit: counts
    path "versions.yml",                emit: versions

    script:
    """
    # Collect R1 (barcodes+UMI) and R2 (cDNA) from all FASTQ directories
    r1=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R1_*.fastq.gz' 2>/dev/null | sort | paste -sd',')
    r2=\$(find ${gex_fastq_dirs.replace(',', ' ')} \\
              -name '${meta.id}*_R2_*.fastq.gz' 2>/dev/null | sort | paste -sd',')

    if [ -z "\$r1" ] || [ -z "\$r2" ]; then
        echo "ERROR: no FASTQs found for ${meta.id} in: ${gex_fastq_dirs}" >&2
        exit 1
    fi

    # Cellbender outputs barcodes with -1 suffix; strip for simpleaf
    sed 's/-1\$//' ${barcodes} > barcodes_clean.txt

    export ALEVIN_FRY_HOME=\${PWD}/.alevin_fry_home
    mkdir -p "\${ALEVIN_FRY_HOME}"
    simpleaf set-paths

    simpleaf quant \\
        --reads1        "\${r1}" \\
        --reads2        "\${r2}" \\
        --threads       ${task.cpus} \\
        --index         ${spliceu_index} \\
        --chemistry     ${simpleaf_chemistry} \\
        --t2g-map       ${spliceu_index}/t2g_3col.tsv \\
        --resolution    cr-like \\
        --unfiltered-pl barcodes_clean.txt \\
        --usa-mode \\
        --output        velocity

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf --version | sed 's/simpleaf //')
    END_VERSIONS
    """
}
```

**Anti-pattern guards**:
- `gex_fastq_dirs` must stay `val(String)` — NOT `path()` (same reason as CELLRANGER_MULTI: NFS paths, multi-flowcell collision avoidance)
- `spliceu_index` is `path` (directory) — Nextflow symlinks entire dir into work dir; `${spliceu_index}/t2g_3col.tsv` resolves correctly inside script
- Do NOT add `--sketch` flag — changes output format
- Do NOT use `metas as ArrayList` in process script — materialise in `.map {}` operators in main.nf

**Verification checklist (Phase 1)**:
- [ ] `nextflow config` parses without syntax errors
- [ ] `velocity/af_quant/` appears in output with `alevin/` subdirectory
- [ ] `versions.yml` present

---

## Phase 2: `get_velocity_chemistry()` Helper

**File**: `main.nf`  
**Location**: Alongside other helper closures (near `detect_sequencer`, `load_si_indexes`).

```groovy
def get_velocity_chemistry(chemistry) {
    switch (chemistry) {
        case 'SC3Pv2':           return '10xv2'
        case 'SC3Pv3':           return '10xv3'
        case 'SC3Pv4':           return '10xv4-3p'
        case 'SC5P':             return '10xv3'
        case 'SC5Pv3':           return '10xv3'
        case 'ARCv1':            return '10xv3'
        default:                 return null   // Flex → skip
    }
}
```

Returns `null` for unsupported chemistries → caller filters these libraries out.

**Verification checklist (Phase 2)**:
- [ ] `get_velocity_chemistry('SC3Pv3')` returns `'10xv3'`
- [ ] `get_velocity_chemistry('Flex-v2-R1')` returns `null`

---

## Phase 3: Wire into `main.nf`

### 3a — Import

Add to GEX module include block:
```groovy
include { SIMPLEAF_VELOCITY } from './modules/qc'
```

### 3b — Build velocity FASTQ channel (before COUNT_GEX)

Extract GEX FASTQs per library, resolve chemistry, key by library_id:

```groovy
// ── RNA Velocity: collect GEX FASTQs ─────────────────────────────────────────
ch_routed.gex
    .filter { meta, fastq_dir ->
        meta.modality == 'GEX' && get_velocity_chemistry(meta.chemistry) != null
    }
    .map { meta, fastq_dir ->
        [ meta.library_id, meta, fastq_dir, get_velocity_chemistry(meta.chemistry) ]
    }
    .groupTuple(by: 0)
    .map { library_id, metas, fastq_dirs, chems ->
        def ml = []; metas.each    { ml << it }
        def dl = []; fastq_dirs.each { dl << it }
        def meta = ml[0] + [library_id: library_id]
        [ library_id, meta, dl.join(','), chems[0] ]
    }
    .set { ch_velocity_fastqs }   // [library_id, meta, fastq_dirs_csv, simpleaf_chemistry]
```

### 3c — Join with cellbender barcodes (after QC_GEX)

After `QC_GEX(ch_gex_outs)`:

```groovy
if (params.run_velocity) {
    // QC_GEX.out.barcodes: [meta, barcodes_csv] — join on library_id
    ch_velocity_fastqs
        .join(
            QC_GEX.out.barcodes.map { meta, barcodes -> [ meta.library_id, barcodes ] },
            by: 0
        )
        .map { library_id, meta, fastq_dirs, chemistry, barcodes ->
            def is_human = meta.species == 'human'
            def idx = file(is_human ? params.spliceu_index_human : params.spliceu_index_mouse)
            [ meta, fastq_dirs, chemistry, barcodes, idx ]
        }
        .set { ch_velocity_ready }

    SIMPLEAF_VELOCITY(
        ch_velocity_ready.map { meta, fqs, chem, bc, idx -> [ meta, fqs, chem, bc ] },
        ch_velocity_ready.map { meta, fqs, chem, bc, idx -> idx }.first()
    )
}
```

**Note**: `spliceu_index` is shared across all samples but species-resolved per-sample. Since `path` inputs must be uniform across all channel items, use `.first()` only if all samples in a run are the same species. If mixed-species runs are supported, switch to per-sample `path` input (remove `.first()`).

**Anti-pattern guards**:
- Do NOT call `get_velocity_chemistry()` inside process `script:` blocks
- Do NOT use `.flatten()` on `metas` from `groupTuple` — use `.each { ml << it }`

### 3d — Add `barcodes` to QC_GEX emit block

**File**: `subworkflows/qc_gex.nf`

```groovy
emit:
    cellbender = CELLBENDER.out.h5
    barcodes   = CELLBENDER.out.barcodes   // ← add this line
    doublets   = SCRUBLET.out.doublets
    vireo      = GENOTYPE.out.donor_ids
    logs       = Channel.empty()
```

**Verification checklist (Phase 3)**:
- [ ] `ch_velocity_fastqs` emits correct tuples (test with `.view()`)
- [ ] Join with barcodes succeeds (no empty channel)
- [ ] No ArrayBag errors in Nextflow log

---

## Phase 4: `nextflow.config` Updates

```groovy
// ── RNA Velocity ──────────────────────────────────────────────────────────────
run_velocity          = true
spliceu_index_human   = "/home/knighto/work/ref/hs/GRCh38-hardmasked-optimised-arc-simpleaf/index"
spliceu_index_mouse   = "/home/knighto/work/ref/mm/GRCm38-hardmasked-optimised-arc-simpleaf/index"
```

Optionally override resources (process_high = 16c/64GB/24h, may be enough):
```groovy
withName: 'SIMPLEAF_VELOCITY' {
    cpus   = 16
    memory = 64.GB
    time   = 12.h
}
```

**Verification checklist (Phase 4)**:
- [ ] `nextflow config -profile slurm` shows both spliceu_index params
- [ ] `--run_velocity false` skips the process

---

## Phase 5: R Import (Seurat)

Output location: `{run_name}_outs/{library_id}/velocity/af_quant/`

```r
library(fishpond)
library(Seurat)

# Load USA mode output
fry <- loadFry(
  fryDir = "path/to/{library_id}/velocity/af_quant",
  outputFormat = "velocity"
)
# fry is a SummarizedExperiment with assays:
#   spliced   — S layer
#   unspliced — U layer
#   ambiguous — A layer

# Extract matrices
spliced   <- assay(fry, "spliced")    # sparse Matrix, genes × cells
unspliced <- assay(fry, "unspliced")

# Add velocity layers to existing Seurat object (cells must match)
seurat_obj[["spliced"]]   <- CreateAssayObject(counts = spliced[, colnames(seurat_obj)])
seurat_obj[["unspliced"]] <- CreateAssayObject(counts = unspliced[, colnames(seurat_obj)])

# Or use scVelo via Python (anndata format)
# import scvelo as scv
# adata = scv.read("velocity/af_quant/")  # requires pyroe/scvelo integration
```

Document this in the pipeline README/docs rather than in CLAUDE.md.

---

## Phase 6: CLAUDE.md Update

Add to:
- **Directory structure**: `{library_id}/velocity/` alongside cellbender/, vireo/
- **Key params**: velocity params block
- **QC section**: SIMPLEAF_VELOCITY description (runs after cellbender; uses cellbender barcodes)
- **Routing rules**: note Flex excluded

---

## Phase 7: End-to-End Verification

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/test.csv \
    --from_fastq  true \
    --fastq_dir   /path/to/test_fastqs \
    --outdir      test_velocity_out \
    --run_name    TEST_VEL
```

- [ ] SIMPLEAF_VELOCITY SLURM job submitted after CELLBENDER completes
- [ ] `{library_id}/velocity/af_quant/` present in output
- [ ] `af_quant/alevin/quants_mat_cols.txt` contains `-S`/`-U`/`-A` gene suffixes
- [ ] Cell count in `quants_mat_rows.txt` matches cellbender barcode count
- [ ] `--run_velocity false` skips process entirely
- [ ] Mouse library routes to `spliceu_index_mouse`
- [ ] fishpond::loadFry() succeeds in R on the output

---

## Known Gaps

1. **Mixed-species runs**: `.first()` on spliceu_index assumes uniform species. If human + mouse in same run, switch to per-sample path input.
2. **`--from-cellranger` mode**: FASTQs unavailable; velocity cannot run. No fix planned — document limitation.
3. **Multi-run libraries** (multiple flowcells): `gex_fastq_dirs` will be comma-separated multiple dirs after groupTuple — the `find` command handles multiple start dirs natively, should work.
4. **Mouse index**: not yet built (Phase 0).
5. **Cellbender `-1` stripping**: if cellbender version changes output format (no `-1` suffix), the `sed` command becomes a no-op (safe).
