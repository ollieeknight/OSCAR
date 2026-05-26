# OSCAR Nextflow Pipeline

Single-cell sequencing pipeline for Romagnani Lab (BIH Charité). Handles scRNA-seq, scATAC-seq, Multiome, DOGMA, ASAP-seq, CITE-seq, Flex (Fixed RNA Profiling), VDJ. Runs on SLURM via Nextflow DSL2 with Apptainer containers.

Legacy bash scripts are in `bash/` for reference only. Nextflow is the primary implementation.

---

## Directory structure

```
.
├── main.nf                    ← entry point; routing, channel logic, helper functions
├── nextflow.config            ← params, profiles (standard / slurm), resource labels
├── assets/
│   ├── example_samplesheet.csv
│   ├── multiqc_config.yaml
│   └── indexes/               ← 10x SI/DI kit CSVs (loaded at runtime for BCL Convert)
│       ├── Single_Index_Kit_GA_Set_A.csv    (SI-GA-*)
│       ├── Single_Index_Kit_NA_Set_A.csv    (SI-NA-*)
│       ├── Dual_Index_Kit_TT_Set_A.csv      (SI-TT-*)
│       ├── Dual_Index_Kit_TN_Set_A.csv      (SI-TN-*)
│       └── truseq_adt_hto.csv               (TruSeq D7xx indices for ADT/HTO)
├── modules/
│   ├── demux.nf      ← GENERATE_SAMPLESHEET, BCLCONVERT, FALCO, MULTIQC
│   ├── count_gex.nf  ← MULTI_CONFIG, CELLRANGER_MULTI
│   ├── count_atac.nf ← CELLRANGER_ATAC
│   ├── count_adt.nf  ← FEATUREMAP, KALLISTO_INDEX, ASAP_TO_KITE,
│   │                    KALLISTO_BUS, BUSTOOLS_CORRECT, BUSTOOLS_SORT, BUSTOOLS_COUNT
│   └── qc.nf         ← CELLBENDER, SCRUBLET, CELLSNP_LITE, VIREO,
│                        AMULET, MGATK2, MACS3
├── subworkflows/
│   ├── demux.nf       ← DEMUX: GENERATE_SAMPLESHEET + BCLCONVERT + FALCO
│   ├── count_gex.nf   ← COUNT_GEX: MULTI_CONFIG → CELLRANGER_MULTI
│   ├── count_atac.nf  ← COUNT_ATAC: thin wrapper around CELLRANGER_ATAC
│   ├── count_adt.nf   ← COUNT_ADT: full ASAP kallisto pipeline
│   ├── qc_gex.nf      ← QC_GEX: CELLBENDER → SCRUBLET → [CELLSNP_LITE → VIREO]
│   └── qc_atac.nf     ← QC_ATAC: AMULET + MGATK2 + MACS3 → [CELLSNP_LITE → VIREO]
├── apptainer/
│   └── mgatk2_docker/ ← Dockerfile + GHA workflow for quay.io/ollieeknight/mgatk2
├── docs/              ← built MkDocs site (GitHub Pages from /docs)
├── bash/              ← legacy bash scripts (reference only)
└── PLAN-*.md          ← implementation plans (safe to delete after execution)
```

---

## Invocation

### Default (BCL → FASTQ → count → QC)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir     /path/to/BCL_folder \
    --outdir      results \
    --run_name    R463
```

### Multi-run (same libraries, multiple flowcells)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet    /path/to/metadata.csv \
    --bcl_dir        /path/to/R462 \
    --extra_bcl_dirs /path/to/R463 \
    --outdir         results \
    --run_name       R462
```

`extra_bcl_dirs` is comma-separated for >2 flowcells. By default the same samplesheet is reused for all BCL dirs; use `--extra_samplesheets` to supply a different samplesheet per extra BCL dir (comma-separated, count must match `extra_bcl_dirs`). FASTQs from each flowcell are demuxed independently then merged by `meta.library_id` before counting.

### From pre-existing FASTQs (`--from_fastq`)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from_fastq  true \
    --fastq_dir   /path/to/fastqs/ \
    --outdir      results
```

FASTQs matched by glob `{fastq_dir}/**/{meta.id}*.fastq.gz`. `meta.id` equals the BCL Convert Sample_ID used during demux.

### From pre-existing cellranger outputs (`--from_cellranger`)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet    /path/to/metadata.csv \
    --from_cellranger true \
    --outs_dir       /path/to/results/ \
    --outdir         results
```

Expects GEX outs at `{outs_dir}/{library_id}/outs` and ATAC outs at `{outs_dir}/{library_id}_ATAC/outs`. Runs QC only.

### Stop after a stage (`--run_until`)

```bash
--run_until FASTQ       # stop after demux; produce FASTQs + falco QC only
--run_until cellranger  # stop after counting; skip QC
```

---

## Samplesheet format

```
assay,experiment_id,historical_number,replicate,modality,chemistry,index_type,index,species,n_donors,adt_file
```

| Column | Values |
|---|---|
| assay | `GEX`, `CITE`, `DOGMA`, `ATAC`, `Multiome`, `ASAP`, `Flex` |
| modality | `GEX`, `ATAC`, `ADT`, `HTO`, `VDJ-T`, `VDJ-B`, `CRISPR`, `GENO` |
| chemistry | `SC3Pv2`, `SC3Pv3`, `SC3Pv4`, `SC5P`, `SC5Pv3`, `ARCv1`, `ATAC`, `Flex-v2-R1`, `Flex-v2-RNA-R2`, `NA` |
| index_type | `SI` or `DI` |
| index | kit code (e.g. `SI-GA-A6`, `SI-TT-A9`) or raw 8-mer sequence |
| species | `human` or `mouse` (case-insensitive; normalised to lowercase) |
| n_donors | integer or `NA` (`NA` → treated as 1; no genotyping) |
| adt_file | stem name of ADT CSV (e.g. `CITE_JIA_adt`); blank for non-ADT rows |

Examples:

```csv
# CITE-seq
CITE,JIA_CD7,4,A,GEX,SC3Pv3,SI,SI-GA-A6,human,3,CITE_JIA_CD7_adt
CITE,JIA_CD7,4,A,ADT,SC3Pv3,SI,AGCTTCAG,human,3,CITE_JIA_CD7_adt

# DOGMA GEX run (GEX + ADT → cellranger multi)
DOGMA,LGL,1,A,GEX,ARCv1,DI,SI-TT-A9,human,2,LGL_ADT
DOGMA,LGL,1,A,ADT,ARCv1,DI,TCTCCGGA,human,2,LGL_ADT

# DOGMA ATAC run (separate invocation)
DOGMA,ILC_IBD,4,A,ATAC,NA,DI,SI-NA-C5,human,4,DOGMA_ILC_IBD_ADT

# Flex (Fixed RNA Profiling)
Flex,IBD,1,A,GEX,Flex-v2-R1,SI,SI-GA-B1,human,1,

# ASAP-seq
ASAP,HDNE,1,A,ATAC,ATAC,DI,SI-NA-C11,human,3,ASAP_HDNE_adt
ASAP,HDNE,1,A,ADT,ATAC,DI,CGAGTAAT,human,3,ASAP_HDNE_adt
```

---

## Meta map

Built in `main.nf:parse_row()`. Every channel tuple starts with a meta map.

```groovy
[
    id:                "${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality}",
    library_id:        "${assay}_${experiment_id}_exp${historical_number}_lib${replicate}",
    assay:             String,   // GEX | CITE | DOGMA | ATAC | Multiome | ASAP | Flex
    experiment_id:     String,
    historical_number: String,
    replicate:         String,
    modality:          String,   // GEX | ATAC | ADT | HTO | VDJ-T | VDJ-B | CRISPR | GENO
    chemistry:         String,
    index_type:        String,   // SI | DI
    index:             String,   // raw value from samplesheet
    index_seqs:        [         // resolved BCL Convert sequences
        is_dual: Boolean,
        rows: [ [i7: String] | [i7: String, i5: String] ]
    ],
    species:           String,   // always lowercase: 'human' | 'mouse'
    n_donors:          Integer,  // NA → 1
    adt_file:          String | null,
    adt_csv_path:      String | null   // absolute path; null when file not found
]
```

`id` = BCL Convert Sample_ID = FASTQ filename prefix.  
`library_id` = cellranger `--id` / output directory name.

### ADT CSV resolution (3-tier, local-first)

1. `{samplesheet_dir}/adt_files/{adt_file}.csv` — co-located with run (preferred)
2. `{samplesheet_dir}/../adt_files/{adt_file}.csv` — parent directory search
3. `{params.adt_files_dir}/{adt_file}.csv` — centralised fallback (requires `--adt_files_dir`)

If none found and `adt_file` is set, pipeline logs a warning and will fail at cellranger multi.

---

## Channel flow

### BCL mode (default)

```
params.samplesheet
    → preflight_samplesheet() → parse_row() → ch_meta

per BCL dir: detect_sequencer(bcl_path) → load_si_indexes()
ch_meta + ch_bcl_dirs → ch_meta_bcl

ch_meta_bcl
    → DEMUX
        → groupTuple by {assay}_{index_type}_{chemistry}_{modality}_{index_len}_{bcl_dir.name}
        → materialise ArrayBag → ArrayList + pre-build is_dual / data_header / data_rows
        → GENERATE_SAMPLESHEET (local executor) → SampleSheet.csv
        → BCLCONVERT → flatMap (match FASTQs back to individual metas by meta.id prefix)
        → ch_fastqs [meta, fastq_dir_string, [matched_fqs]]
        → FALCO → ch_falco_reports

ch_fastqs
    → .branch {
        gex:      modality ∈ {GEX,ADT,HTO,VDJ-T,VDJ-B,CRISPR} && assay != ASAP
                  OR assay == Flex
        atac:     modality == ATAC
        asap_adt: assay == ASAP && modality ∈ {ADT,HTO}
        skip:     GENO / else
      }

ch_routed.gex
    → map [library_id, meta, fastq_dir]
    → groupTuple(by: 0)
    → map: materialise ArrayBag, deduplicate metas by modality, collect unique dirs, resolve adt_csv
    → COUNT_GEX:
        → MULTI_CONFIG (local executor): build multi_config.csv + validate read counts
        → CELLRANGER_MULTI: cellranger multi --csv multi_config.csv
    → ch_gex_outs [library_id, metas, outs_dir]
    → QC_GEX:
        → CELLBENDER (GPU)
        → (if n_donors > 1 && species == human) CELLSNP_LITE → VIREO

ch_routed.atac
    → map [library_id, meta, fastq_dir] → groupTuple(by: 0)
    → COUNT_ATAC → CELLRANGER_ATAC
    → ch_atac_outs [meta, outs_dir]
    → QC_ATAC:
        → AMULET  ─┐
        → MGATK2   ├─ (all three always run, in parallel)
        → MACS3   ─┘
        → (if n_donors > 1 && species == human) CELLSNP_LITE → VIREO

ch_atac_outs (ASAP only) + ch_routed.asap_adt
    → join by library_id
    → COUNT_ADT:
        → FEATUREMAP → KALLISTO_INDEX
        → ASAP_TO_KITE → KALLISTO_BUS → BUSTOOLS_CORRECT → BUSTOOLS_SORT → BUSTOOLS_COUNT

ch_falco_reports + QC_GEX.out.logs → MULTIQC
```

---

## Demux

Container: `quay.io/nf-core/bclconvert:4.4.6`

### Sequencer auto-detection

`detect_sequencer(bcl_path)` reads `<Instrument>` from `RunInfo.xml` per BCL dir:

| Prefix | Instrument | i5 orientation |
|---|---|---|
| `VH` | NovaSeq X / X Plus | forward → `novaseq_x` |
| `NDX` | NextSeq 2000/1000 | forward → `novaseq_x` |
| `FS` | iSeq 100 | forward → `novaseq_x` |
| `A`, `LH` | NovaSeq 6000 | reverse-complement → `novaseq6000` |
| `NB`, `NS`, `MN` | NextSeq 550/500, MiniSeq | reverse-complement → `novaseq6000` |

Falls back to `params.sequencer` (default: `novaseq_x`) if RunInfo.xml absent or instrument unrecognised. Each BCL dir in a multi-run is detected independently — a pipeline run can mix instruments.

### Index resolution

SI/DI kit codes resolved at startup by `load_si_indexes()` from `assets/indexes/` CSVs:

- `SI-GA-*` / `SI-NA-*` → 4 single-index rows (8 bp each), `is_dual = false`
- `SI-TT-*` / `SI-TN-*` → 1 dual-index row (i7 + i5), `is_dual = true`
- Raw 8-mer → 1 single-index row, `is_dual = false`

**SI-on-DI correction**: When `index_type=DI` but the actual index is single-index, `get_override_cycles()` corrects: I1 → `I{len}N*`, I2 → `N*` (fully masked).

### OverrideCycles

`get_base_mask()` defines the mask table; `get_override_cycles()` converts to BCL Convert format (`;` separator, uppercase `N`) and resolves `*` wildcards to exact cycle counts at runtime by reading `NumCycles` from `RunInfo.xml`.

### Demux key

`{assay}_{index_type}_{chemistry}_{modality}_{index_len}_{bcl_dir.name}` — the BCL dir name suffix prevents work-dir collision when the same library group is demuxed across multiple flowcells.

### GENERATE_SAMPLESHEET

Runs on `local` executor (not SLURM) because it only writes a CSV — no compute needed. Receives `is_dual`, `data_header`, and `data_rows` as pre-built strings from the subworkflow `.map {}` to avoid Groovy collection ops inside a process script block.

---

## Counting

### GEX / CITE / DOGMA-GEX / Multiome-GEX / Flex → COUNT_GEX

Container: `quay.io/nf-core/cellranger:10.0.0`

COUNT_GEX is a two-step subworkflow:

**Step 1 — MULTI_CONFIG** (local executor): generates `multi_config.csv`.
- Config sections built in the subworkflow `.map {}` as a `List<String>`, joined and passed as a pre-built string — no Groovy collection ops in the process script block.
- Per (fastq_dir × modality): reads first 40000 lines of R1 via `zcat | head | awk` to count reads. Pairs with <10000 reads are excluded — threshold matches cellranger's own chemistry auto-detection minimum (TXRNGR10001). Flowcells contributing fewer reads would fail cellranger regardless.
- Threshold: `min_reads = 10000` (Groovy variable in subworkflow, easy to tune).

**Step 2 — CELLRANGER_MULTI**: pure bash, no Groovy logic. Receives `path(multi_config)` and runs `cellranger multi --csv multi_config`.

**`[gene-expression]` section chemistry:**

| Assay | Added line |
|---|---|
| DOGMA, Multiome | `chemistry,ARC-v1` |
| Flex | `chemistry,{meta.chemistry}` (e.g. `Flex-v2-R1`) |
| Everything else | _(no chemistry line)_ |

ADT CSV: referenced in `[feature]` section. Uses `NO_FILE` sentinel when `--adt_files_dir` not set or file not found — skips `[feature]` section entirely.

Output: `{outdir}/{run_name}_outs/{library_id}/outs`

### ATAC / DOGMA-ATAC / Multiome-ATAC / ASAP-ATAC → COUNT_ATAC

Container: `quay.io/nf-core/cellranger-atac:2.1.0`

DOGMA/Multiome ATAC adds `--chemistry ARC-v1`.

Output: `{outdir}/{run_name}_outs/{library_id}_ATAC/outs`

### ASAP ADT → COUNT_ADT (kallisto pipeline)

Triggered by channel join with CELLRANGER_ATAC output — runs only after ATAC counting completes for the matching library.

```
FEATUREMAP    (asap_to_kite)  : kite featuremap → FeaturesMismatch.fa + .t2g
KALLISTO_INDEX (kallisto)     : kallisto index -k 15
ASAP_TO_KITE   (asap_to_kite) : converts ATAC barcodes → GEX barcodes
KALLISTO_BUS   (kallisto)     : -x 0,0,16:0,16,26:1,0,0
BUSTOOLS_CORRECT (bustools)   : correct against params.atac_whitelist
BUSTOOLS_SORT
BUSTOOLS_COUNT                : --genecounts
```

Output: `{library_id}_ATAC/ADT/`

---

## QC

### GEX QC — QC_GEX

Input: `[library_id, metas, outs_dir]` from COUNT_GEX.

```
CELLBENDER  (GPU)   — ambient RNA removal
  ↓
SCRUBLET            — GEX doublet detection (outputs doublets.csv metadata)
  ↓ (if n_donors > 1 && species == human)
CELLSNP_LITE  : pileup from per_sample_outs BAM (mode='gex', UMI-tagged)
VIREO         : probabilistic donor demultiplexing
```

CELLBENDER uses the raw feature-barcode matrix (`raw_feature_bc_matrix.h5`). Requires GPU node.

Output channels: `cellbender` (h5), `doublets` (doublets.csv), `vireo` (donor_ids.tsv, empty if single-donor).

### ATAC QC — QC_ATAC

Input: `[meta, outs_dir]` from COUNT_ATAC.

```
AMULET  ─┐
MGATK2   ├─ all three run always, in parallel
MACS3   ─┘
  ↓ (if n_donors > 1 && species == human)
CELLSNP_LITE  : pileup from possorted_bam.bam (mode='atac', no UMI tag)
VIREO         : probabilistic donor demultiplexing
```

**MACS3** settings follow ENCODE scATAC recommendations: `--nomodel --shift -75 --extsize 150 --keep-dup all --nolambda`. FDR threshold: `params.macs3_qvalue` (default `0.05`). Streams fragments via process substitution to avoid large temp files.  
**AMULET** uses species-matched autosome and restriction/repeat lists from `/opt/AMULET/`.  
**MGATK2** runs mitochondrial genotyping on `possorted_bam.bam`.

Output channels: `amulet` (summary), `mgatk` (results dir), `peaks` (MACS3 peaks dir), `vireo` (donor_ids.tsv).

---

## Output structure

```
results/
├── {run_name}_fastq/                  ← BCL Convert FASTQs (BCL mode only)
│   └── falco/                         ← per-FASTQ quality reports
│
├── {run_name}_outs/
│   ├── {library_id}/                  ← GEX / CITE / DOGMA-GEX / Flex output
│   │   └── outs/
│   │       ├── per_sample_outs/
│   │       ├── multi/
│   │       ├── cellbender/            ← output.h5, output_cell_barcodes.csv
│   │       ├── scrublet/              ← doublets.csv (GEX doublet scores)
│   │       └── vireo/                 ← donor_ids.tsv, cellsnp/ (multi-donor only)
│   │
│   └── {library_id}_ATAC/             ← ATAC / DOGMA-ATAC output
│       └── outs/
│           ├── ADT/                   ← ASAP kallisto counts (ASAP only)
│           ├── AMULET/                ← MultipletSummary.txt, MultipletBarcodes.txt
│           ├── mgatk2/                ← mgatk2_out/
│           ├── peaks/                 ← MACS3 peak files
│           └── vireo/                 ← donor_ids.tsv, cellsnp/ (multi-donor only)
│
└── multiqc/
    └── multiqc_report.html
```

---

## Containers

```groovy
// 10x Genomics (nf-core public images)
container_bclconvert      = "quay.io/nf-core/bclconvert:4.4.6"
container_cellranger      = "quay.io/nf-core/cellranger:10.0.0"
container_cellranger_atac = "quay.io/nf-core/cellranger-atac:2.1.0"

// ATAC QC (custom/lab images)
container_amulet = "quay.io/cellgeni/amulet:1.1"
container_mgatk  = "quay.io/ollieeknight/mgatk2:latest"   // built from apptainer/mgatk2_docker/

// ASAP/ADT pipeline (public lab image)
container_asap = "quay.io/ollieeknight/asap_to_kite:main"

// GEX QC
container_cellbender = "quay.io/biocontainers/cellbender:0.3.2--pyhdfd78af_0"

// ATAC peak calling
container_macs3 = "quay.io/biocontainers/macs3:3.0.1--py311h320fe9a_1"

// Genotyping / demultiplexing
container_cellsnp = "quay.io/biocontainers/cellsnp-lite:1.2.3--ha0c3a46_6"
container_vireo   = "quay.io/biocontainers/vireosnp:0.5.9--pyh7e72e81_0"

// QC / reporting
container_falco   = "quay.io/biocontainers/falco:1.3.0--h3be2455_0"
container_multiqc = "quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0"

// ASAP kallisto pipeline
container_kallisto = "quay.io/biocontainers/kallisto:0.51.1--h7d033e6_0"
container_bustools = "quay.io/biocontainers/bustools:0.45.1--h9f5acd7_0"
```

---

## Key params

```groovy
// ── Input ──────────────────────────────────────────────────────────────────────
samplesheet     = null     // path to metadata.csv
bcl_dir         = null     // path to BCL folder (RunInfo.xml inside)
outdir          = 'results'
run_name        = null     // prefix for _fastq / _outs dirs; auto-derived from bcl_dir (R463_bcl → R463) if not set

// Entry point (exactly one active at a time)
from_fastq      = false    // start from pre-existing FASTQs
from_cellranger = false    // start from cellranger outs (QC only)
fastq_dir       = null     // required with from_fastq
outs_dir        = null     // required with from_cellranger

// Pipeline control
run_until       = null     // 'FASTQ' | 'cellranger' | null (full run)

// Multi-run
extra_bcl_dirs      = null // comma-separated additional BCL dirs
extra_samplesheets  = null // comma-separated samplesheets aligned with extra_bcl_dirs

// ADT files
adt_files_dir   = null     // centralised fallback dir for {adt_file}.csv
                            // (local search relative to samplesheet is tried first)

// ── Sequencer ─────────────────────────────────────────────────────────────────
sequencer       = 'novaseq_x'  // fallback when RunInfo.xml absent or unrecognised
                                // auto-detected per BCL dir from RunInfo.xml in BCL mode

// ── References ────────────────────────────────────────────────────────────────
ref_human     = "/sc-projects/.../GRCh38-hardmasked-optimised-arc"
ref_vdj_human = "/sc-projects/.../GRCh38-IMGT-VDJ-2024"
ref_mouse     = "/sc-projects/.../GRCm38-hardmasked-optimised-arc"
ref_vdj_mouse = "/sc-projects/.../GRCm38-IMGT-VDJ-2024"
snp_vcf       = "/sc-projects/.../genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
atac_whitelist = "/home/knighto/bin/cellranger-atac-2.2.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz"

// ── QC tuning ─────────────────────────────────────────────────────────────────
macs3_qvalue           = 0.05   // FDR threshold for peak calling
```

---

## SLURM resource labels

```
process_low     →  2c /  8GB /  4h
process_medium  →  8c / 32GB / 12h
process_high    → 16c / 64GB / 24h
process_gpu     → 16c / 96GB / 24h  + queue=gpu, --gres=gpu:1, --nv (Apptainer --nv)

withName: GENERATE_SAMPLESHEET|MULTI_CONFIG → local executor (no SLURM job; CSV generation only)
withName: BCLCONVERT                      → 16c / 64GB / 24h
withName: CELLRANGER_MULTI|CELLRANGER_ATAC  → 64c / 128GB / 48h
withName: CELLSNP_LITE|VIREO                → 32c /  64GB / 48h
withName: MGATK2                            → 32c / 128GB / 48h
withName: KALLISTO_BUS                      → 16c /  96GB / 48h
```

Process labels in modules:
- `AMULET`, `MGATK2`, `MACS3` → `process_medium` (MGATK2 overridden above)
- `CELLBENDER` → `process_gpu`
- `CELLRANGER_MULTI`, `CELLRANGER_ATAC` → `process_high` (overridden above)

---

## Routing rules

| Assay | Modality | Route |
|---|---|---|
| GEX, CITE | GEX, ADT, HTO, VDJ-T, VDJ-B, CRISPR | CELLRANGER_MULTI → QC_GEX |
| DOGMA | GEX, ADT, HTO | CELLRANGER_MULTI → QC_GEX |
| DOGMA | ATAC | CELLRANGER_ATAC → QC_ATAC |
| Multiome | GEX, ADT | CELLRANGER_MULTI → QC_GEX |
| Multiome | ATAC | CELLRANGER_ATAC → QC_ATAC |
| Flex | GEX | CELLRANGER_MULTI → QC_GEX |
| ATAC | ATAC | CELLRANGER_ATAC → QC_ATAC |
| ASAP | ATAC | CELLRANGER_ATAC → QC_ATAC + COUNT_ADT |
| ASAP | ADT, HTO | COUNT_ADT (kallisto) |
| Any | GENO | skipped |

**DOGMA ADT never goes to kallisto.** It is handled by cellranger multi on the GEX run.  
**DOGMA GEX + ATAC are separate invocations.** Integration done downstream in R, not in Nextflow.

---

## Design constraints

### Groovy collection ops: map blocks vs process script blocks

`groupTuple` emits `nextflow.util.ArrayBag` — Groovy's internal grouped-list type. Calling `.any {}`, `.collect {}`, `.collectMany {}`, `.flatten()`, or `as ArrayList` on an ArrayBag inside a process `script:` block throws `UnsupportedOperationException`.

**Rule**: All collection logic (iterating metas, building strings, pre-computing values) must run inside workflow `.map {}` operators — never inside process `script:` blocks. Process script blocks receive only plain types: `String`, `Boolean`, `Integer`, pre-joined strings.

Pattern for safe materialisation in a `.map {}`:
```groovy
.map { key, metas, ... ->
    def ml = []
    metas.each { m -> ml << m }   // ← safe: iterates ArrayBag, appends to ArrayList
    // now ml is java.util.ArrayList — safe for .any {}, .collect {}, .collectMany {}
}
```

This pattern is used in:
- `subworkflows/count_gex.nf` — materialises metas before MULTI_CONFIG
- `subworkflows/demux.nf` — materialises metas and pre-builds `is_dual`, `data_header`, `data_rows` before GENERATE_SAMPLESHEET

### val(fastq_dirs) — directory strings, not staged files

GEX counting passes fastq directories as `val(String)` rather than staging files via `path()`. This prevents filename collisions when the same library is sequenced on multiple flowcells (cellranger produces identical `*_S1_*.fastq.gz` names). Cellranger reads from NFS paths directly.

### FASTQ placeholder filtering

BCL Convert produces near-empty placeholder FASTQs for samples absent from a flowcell. These can exceed 10 MB but contain <100 real reads. MULTI_CONFIG filters them by read count (threshold: 1000 reads) rather than file size:

```bash
r1=$(find "${dir}" -maxdepth 2 -name "${m.id}*_R1_*.fastq.gz" -print -quit)
n_reads=0; [ -n "$r1" ] && n_reads=$(zcat "$r1" | head -n 10000 | awk 'NR%4==1' | wc -l)
[ "$n_reads" -ge 2500 ] && echo "${m.id},${dir},feature_type" >> multi_config.csv || true
```

### No Nextflow staging for ATAC FASTQs

`CELLRANGER_ATAC` also receives `val(fastq_dirs)` and passes `--fastqs` as directory strings. Same reason as GEX.

---

## What NOT to do

- Do not add columns to the samplesheet — format must match `bash/metadata.csv`
- Do not merge DOGMA GEX + ATAC outputs in Nextflow — R integration downstream
- Do not route DOGMA/Flex ADT to kallisto — goes to cellranger multi
- Do not put Groovy collection logic (`.any {}`, `.collect {}`, `.flatten()`, `as ArrayList`) inside process `script:` blocks — use `.map {}` in the subworkflow instead
- Do not use `[metas].flatten()` or `metas as ArrayList` on a value that came from `groupTuple` — always iterate with `.each { ml << it }` into a `[]` literal
- Do not invent Nextflow operators — confirmed: `.branch`, `.join`, `.groupTuple`, `.combine`, `.flatten`, `.multiMap`, `.mix`, `.collect`, `.map`, `.filter`, `.set`, `groupKey()`
- Every process must emit `versions.yml`
- `n_donors=NA` → treat as 1; skip genotyping
- `adt_csv_path` is null when no CSV found — CELLRANGER_MULTI receives `file('NO_FILE')` sentinel and omits `[feature]` section
- Do not hardcode `params.sequencer` — it auto-detects from RunInfo.xml; `params.sequencer` is only the fallback
