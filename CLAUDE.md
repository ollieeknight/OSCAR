# OSCAR Nextflow Pipeline

Single-cell sequencing pipeline for Romagnani Lab (BIH Charité). Handles scRNA-seq, scATAC-seq, Multiome, DOGMA, ASAP-seq, CITE-seq, VDJ. Runs on SLURM via Nextflow DSL2 with Apptainer containers.

Source bash scripts live in `../bash/` (01–05 + functions.sh). This Nextflow port is complete and replaces them.

---

## Directory structure

```
nextflow/
├── main.nf                    ← entry point; all routing and channel logic
├── nextflow.config            ← params, profiles (standard / slurm), resource labels
├── assets/
│   ├── example_samplesheet.csv
│   ├── multiqc_config.yaml
│   └── indexes/               ← 10x SI/DI kit CSV files (loaded at runtime for BCL Convert)
│       ├── Single_Index_Kit_GA_Set_A.csv    (SI-GA-*)
│       ├── Single_Index_Kit_NA_Set_A.csv    (SI-NA-*)
│       ├── Dual_Index_Kit_TT_Set_A.csv      (SI-TT-*)
│       ├── Dual_Index_Kit_TN_Set_A.csv      (SI-TN-*)
│       └── truseq_adt_hto.csv               (TruSeq D7xx indices for ADT/HTO)
├── modules/
│   ├── demux.nf      ← BCL_TO_FASTQ, FALCO, MULTIQC
│   ├── count_gex.nf  ← CELLRANGER_MULTI
│   ├── count_atac.nf ← CELLRANGER_ATAC
│   ├── count_adt.nf  ← FEATUREMAP, KALLISTO_INDEX, ASAP_TO_KITE,
│   │                    KALLISTO_BUS, BUSTOOLS_CORRECT, BUSTOOLS_SORT, BUSTOOLS_COUNT
│   └── qc.nf         ← CELLBENDER, CELLSNP_LITE, VIREO, AMULET, MGATK2
├── subworkflows/
│   ├── demux.nf       ← DEMUX: BCL_TO_FASTQ + FALCO
│   ├── count_gex.nf   ← COUNT_GEX: thin wrapper around CELLRANGER_MULTI
│   ├── count_atac.nf  ← COUNT_ATAC: thin wrapper around CELLRANGER_ATAC
│   ├── count_adt.nf   ← COUNT_ADT: full ASAP kallisto pipeline
│   ├── qc_gex.nf      ← QC_GEX: CELLBENDER → CELLSNP_LITE → VIREO
│   └── qc_atac.nf     ← QC_ATAC: AMULET + MGATK2 → CELLSNP_LITE → VIREO
└── apptainer/
    └── mgatk2_docker/ ← Dockerfile + GHA workflow for quay.io/ollieeknight/mgatk2
```

---

## Invocation

### Default (BCL → FASTQ → count → QC)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --bcl_dir     /path/to/BCL_folder \
    --outdir      results \
    --adt_files_dir /path/to/adt_csvs/
```

### Multi-run (same libraries, two flowcells)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet   /path/to/metadata.csv \
    --bcl_dir       /path/to/R462 \
    --extra_bcl_dirs /path/to/R463 \
    --outdir        results \
    --adt_files_dir /path/to/adt_csvs/
```

`extra_bcl_dirs` is comma-separated for >2 flowcells. Same samplesheet is used for all BCL dirs. FASTQs from each flowcell are demuxed independently then merged by `meta.library_id` before counting.

### From pre-existing FASTQs (`--from-fastq`)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet /path/to/metadata.csv \
    --from-fastq  true \
    --fastq_dir   /path/to/fastqs/ \
    --outdir      results \
    --adt_files_dir /path/to/adt_csvs/
```

FASTQs matched by glob `{fastq_dir}/**/{meta.id}*.fastq.gz`. `meta.id` equals the BCL Convert Sample_ID used during demux.

### From pre-existing cellranger outputs (`--from-cellranger`)

```bash
nextflow run main.nf -profile slurm \
    --samplesheet    /path/to/metadata.csv \
    --from-cellranger true \
    --outs_dir       /path/to/results/ \
    --outdir         results
```

Expects GEX outs at `{outs_dir}/{library_id}/outs` and ATAC outs at `{outs_dir}/{library_id}_ATAC/outs`. Runs QC only; skips demux and counting.

---

## Samplesheet format

Identical to current OSCAR `metadata.csv` — no new columns.

```
assay,experiment_id,historical_number,replicate,modality,chemistry,index_type,index,species,n_donors,adt_file
```

| Column | Values / notes |
|---|---|
| assay | GEX, CITE, DOGMA, ATAC, Multiome, ASAP |
| modality | GEX, ATAC, ADT, HTO, VDJ-T, VDJ-B, CRISPR, GENO |
| chemistry | SC3Pv2/v3/v4, SC5P, SC5Pv3, ARCv1, ATAC, NA |
| index_type | SI or DI |
| index | SI-GA-A6, SI-NA-C5, SI-TT-A9, or raw 8-mer sequence |
| species | human or mouse (case-insensitive; normalised to lowercase internally) |
| n_donors | integer or NA (NA treated as 1; no genotyping run) |
| adt_file | stem name of ADT CSV in `--adt_files_dir` (e.g. `CITE_JIA_adt`); null for non-ADT rows |

Real examples:

```csv
# CITE-seq: GEX + ADT on GEX-type run
CITE,JIA_CD7,4,A,GEX,SC3Pv3,SI,SI-GA-A6,human,3,CITE_JIA_CD7_adt
CITE,JIA_CD7,4,A,ADT,SC3Pv3,SI,AGCTTCAG,human,3,CITE_JIA_CD7_adt

# DOGMA GEX run (GEX + ADT → cellranger multi, same as CITE)
DOGMA,LGL,1,A,GEX,ARCv1,DI,SI-TT-A9,human,2,LGL_ADT
DOGMA,LGL,1,A,ADT,ARCv1,DI,TCTCCGGA,human,2,LGL_ADT

# DOGMA ATAC run (separate invocation → cellranger-atac only)
DOGMA,ILC_IBD,4,A,ATAC,NA,DI,SI-NA-C5,human,4,DOGMA_ILC_IBD_ADT

# ASAP-seq: ATAC → cellranger-atac, ADT → kallisto pipeline
ASAP,HDNE,1,A,ATAC,ATAC,DI,SI-NA-C11,Human,3,ASAP_HDNE_adt
ASAP,HDNE,1,A,ADT,ATAC,DI,CGAGTAAT,Human,3,ASAP_HDNE_adt
```

---

## Meta map

Built in `main.nf:parse_row()`. Every channel tuple starts with a meta map.

```groovy
[
    id:               "${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality}",
    library_id:       "${assay}_${experiment_id}_exp${historical_number}_lib${replicate}",
    assay:            String,          // GEX | CITE | DOGMA | ATAC | Multiome | ASAP
    experiment_id:    String,
    historical_number: String,
    replicate:        String,
    modality:         String,          // GEX | ATAC | ADT | HTO | VDJ-T | VDJ-B | CRISPR | GENO
    chemistry:        String,
    index_type:       String,          // SI | DI
    index:            String,          // raw value from samplesheet
    index_seqs:       [                // resolved BCL Convert sequences
        is_dual: Boolean,
        rows: [ [i7: String] | [i7: String, i5: String] ]
    ],
    species:          String,          // always lowercase: 'human' | 'mouse'
    n_donors:         Integer,         // NA → 1
    adt_file:         String | null,   // stem name only
    adt_csv_path:     String | null    // absolute path: {adt_files_dir}/{adt_file}.csv
]
```

`id` = BCL Convert Sample_ID = FASTQ filename prefix. `library_id` = cellranger `--id` / output directory name.

---

## Channel flow (BCL mode)

```
params.samplesheet
    → splitCsv → parse_row() → ch_meta

ch_meta + ch_bcl_dirs
    → DEMUX
        → groupTuple by {assay}_{index_type}_{chemistry}_{modality}
        → combine(ch_bcl_dirs)                        [cross-product for multi-run]
        → BCL_TO_FASTQ                                [one job per group × flowcell]
        → flatMap (match FASTQs back to individual metas by meta.id in filename)
        → ch_fastqs [meta, [fastqs]]
        → FALCO → collect → ch_falco_reports

ch_fastqs
    → .branch {
        gex:      modality ∈ {GEX,ADT,HTO,VDJ-T,VDJ-B,CRISPR} && assay != ASAP
        atac:     modality == ATAC
        asap_adt: assay == ASAP && modality ∈ {ADT,HTO}
        skip:     GENO / else
      }

ch_routed.gex
    → map [library_id, meta, [fastqs]]
    → groupTuple(by: 0) + flatten fqs + resolve adt_csv
    → COUNT_GEX → CELLRANGER_MULTI
    → ch_gex_outs [library_id, metas, outs_dir]
    → QC_GEX → CELLBENDER (GPU) → [CELLSNP_LITE → VIREO] if n_donors > 1

ch_routed.atac
    → COUNT_ATAC → CELLRANGER_ATAC
    → ch_atac_outs [meta, outs_dir]
    → QC_ATAC → AMULET + MGATK2 → [CELLSNP_LITE → VIREO] if n_donors > 1

ch_atac_outs (ASAP only) + ch_routed.asap_adt
    → join by library_id
    → COUNT_ADT → FEATUREMAP → KALLISTO_INDEX
                              → ASAP_TO_KITE → KALLISTO_BUS
                                             → BUSTOOLS_CORRECT → SORT → COUNT

DEMUX.out.falco_reports + QC_GEX.out.logs
    → MULTIQC
```

---

## Demux: BCL Convert V2

Uses `quay.io/nf-core/bclconvert:4.4.6_1`. Generates V2 SampleSheet at runtime inside `BCL_TO_FASTQ`.

### Index resolution

SI/DI kit codes (e.g. `SI-GA-A6`, `SI-TT-A9`) are resolved to actual sequences at pipeline startup by `load_si_indexes()` reading four CSV files from `assets/indexes/`. Direct 8-mer sequences pass through as-is. Resolution stored in `meta.index_seqs`.

- SI-GA/SI-NA kits → 4 single-index rows (each 8 bp), `is_dual = false`
- SI-TT/SI-TN kits → 1 dual-index row (i7 + i5), `is_dual = true`
- Raw 8-mer → 1 single-index row, `is_dual = false`
- i5 orientation: `novaseq6000` (default) uses reverse-complement column; `novaseq_x` uses forward column (set `params.sequencer`)

### OverrideCycles

Base mask table ported from `functions.sh` into `get_base_mask()` in `modules/demux.nf`. Converted from bcl2fastq format (`,`, lowercase `n`) to BCL Convert format (`;`, uppercase `N`) by `get_override_cycles()`.

**SI-on-DI correction**: When a library has `index_type=DI` (4-read flow cell) but the actual index is single-index (SI-GA kit or raw 8-mer), the base mask has `I10N*;I10N*` for both index positions. `get_override_cycles` detects `!index_seqs.is_dual` and corrects:
- I1 position → `I{actual_seq_len}N*` (clamped to 8 bp for SI-GA)
- I2 position → `N*` (fully masked, no i5 index present)
- Y positions (data reads and ATAC cell barcodes) are never touched

### Demux key and multi-run

Key = `{assay}_{index_type}_{chemistry}_{modality}_{bcl_dir.name}`. The BCL dir name suffix prevents work-dir collision when the same library group is demuxed from multiple flowcells. After demux, FASTQs from all flowcells flow into the same channel and are merged naturally by `groupTuple(by: library_id)` before counting.

---

## Counting

### GEX / CITE / DOGMA-GEX / Multiome-GEX → CELLRANGER_MULTI

Container: `quay.io/nf-core/cellranger:10.0.0`

All modalities for the same `library_id` are grouped together. FASTQs are staged flat in the work dir. The cellranger multi config uses `fastq_id = meta.id` and `fastqs = .` — cellranger finds each modality's FASTQs by filename prefix matching.

ADT CSV: staged as `path(adt_csv)` input when `params.adt_files_dir` is set; skipped (`NO_FILE` sentinel) otherwise. Referenced directly in the `[feature]` section.

DOGMA and Multiome add `chemistry,ARC-v1` to the `[gene-expression]` section.

Output: `{outdir}/{library_id}/outs`

### ATAC / DOGMA-ATAC / Multiome-ATAC / ASAP-ATAC → CELLRANGER_ATAC

Container: `quay.io/nf-core/cellranger-atac:2.1.0`

DOGMA ATAC adds `--chemistry ARC-v1`.

Output: `{outdir}/{library_id}_ATAC/outs`

### ASAP ADT → kallisto pipeline (COUNT_ADT)

Triggered by channel dependency: runs only after `CELLRANGER_ATAC` completes for the matching ASAP library. Not SLURM-dependency based — purely Nextflow channel ordering.

```
FEATUREMAP    (container_asap)         : kite featuremap → FeaturesMismatch.fa + .t2g
KALLISTO_INDEX (container_kallisto)    : kallisto index -k 15
ASAP_TO_KITE   (container_asap)        : converts ATAC barcodes → GEX barcodes
KALLISTO_BUS   (container_kallisto)    : -x 0,0,16:0,16,26:1,0,0
BUSTOOLS_CORRECT (container_bustools)  : correct against params.atac_whitelist
BUSTOOLS_SORT
BUSTOOLS_COUNT                         : --genecounts → {library_id}_ATAC/ADT/
```

ADT CSV resolved from `adt_meta.adt_csv_path` (set in `parse_row` via `--adt_files_dir`).

---

## QC

### GEX QC (QC_GEX)

```
CELLBENDER (process_gpu)       : ambient RNA removal
  ↓ (if n_donors > 1)
CELLSNP_LITE                   : genotype cells from BAM
VIREO                          : donor demultiplexing
```

CellBender container: `params.container_cellbender` (custom .sif, GPU required).

### ATAC QC (QC_ATAC)

```
AMULET (process_high)          : doublet detection
MGATK2 (process_high)          : mitochondrial genotyping
  ↓ (if n_donors > 1)
CELLSNP_LITE → VIREO
```

AMULET container: `quay.io/cellgeni/amulet:1.1`
mgatk2 container: `quay.io/ollieeknight/mgatk2:latest` (custom, built from `apptainer/mgatk2_docker/`)

---

## Output structure

```
results/
├── {library_id}/                          ← GEX/CITE/DOGMA-GEX output
│   └── outs/
│       ├── per_sample_outs/
│       └── multi/
│
├── {library_id}_ATAC/                     ← ATAC/DOGMA-ATAC output
│   └── outs/
│       ├── ADT/                           ← ASAP kallisto counts (ASAP only)
│
└── multiqc/
    └── multiqc_report.html
```

QC outputs are written inside the library outs dir by each QC process (`publishDir` within each process).

---

## Containers

```groovy
// 10x tools (nf-core public images)
container_bclconvert      = "quay.io/nf-core/bclconvert:4.4.6_1"
container_cellranger      = "quay.io/nf-core/cellranger:10.0.0"
container_cellranger_atac = "quay.io/nf-core/cellranger-atac:2.1.0"

// QC tools
container_amulet   = "quay.io/cellgeni/amulet:1.1"
container_mgatk    = "quay.io/ollieeknight/mgatk2:latest"

// Custom .sif (no public equivalent)
container_asap       = "/home/knighto/bin/apptainer/cache/oscar-asap.sif"
container_cellbender = "/home/knighto/bin/apptainer/cache/oscar-cellbender.sif"

// Public biocontainers (pinned)
container_falco    = "quay.io/biocontainers/falco:1.3.0--h3be2455_0"
container_multiqc  = "quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0"
container_kallisto = "quay.io/biocontainers/kallisto:0.51.1--h7d033e6_0"
container_bustools = "quay.io/biocontainers/bustools:0.45.1--h9f5acd7_0"
container_cellsnp  = "quay.io/biocontainers/cellsnp-lite:1.2.3--ha0c3a46_6"
container_vireo    = "quay.io/biocontainers/vireosnp:0.5.9--pyh7e72e81_0"
```

---

## Key params

```groovy
// Input
samplesheet     = null     // path to metadata.csv
bcl_dir         = null     // path to BCL folder (RunInfo.xml inside)
outdir          = 'results'

// Entry point (default is BCL; exactly one of the below may be set)
from_fastq      = false    // skip demux; load FASTQs from fastq_dir
from_cellranger = false    // skip demux + counting; run QC only from outs_dir
fastq_dir       = null     // required when from_fastq = true
outs_dir        = null     // required when from_cellranger = true

// Multi-run
extra_bcl_dirs  = null     // comma-separated additional BCL dirs

// ADT files
adt_files_dir   = null     // directory containing {adt_file}.csv files

// References (human)
ref_human       = "/sc-projects/.../GRCh38-hardmasked-optimised-arc"
ref_vdj_human   = "/sc-projects/.../GRCh38-IMGT-VDJ-2024"

// References (mouse)
ref_mouse       = "/sc-projects/.../GRCm38-hardmasked-optimised-arc"
ref_vdj_mouse   = "/sc-projects/.../GRCm38-IMGT-VDJ-2024"

// QC
snp_vcf         = "/sc-projects/.../genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
atac_whitelist  = "/home/knighto/bin/cellranger-atac-2.2.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz"

// Sequencer (controls i5 orientation for dual-index kits)
sequencer       = 'novaseq6000'   // or 'novaseq_x'
```

---

## SLURM resource labels

```
process_low     →  2c /  8GB /  4h
process_medium  →  8c / 32GB / 12h
process_high    → 16c / 64GB / 24h
process_gpu     → 16c / 96GB / 24h  + queue=gpu, --gres=gpu:1, --nv

withName: CELLRANGER_MULTI|CELLRANGER_ATAC → 64c / 128GB / 48h
withName: BCL_TO_FASTQ                     → 16c /  32GB / 12h
withName: CELLSNP_LITE|VIREO               → 32c /  64GB / 48h
withName: MGATK2                           → 32c / 128GB / 48h
withName: KALLISTO_BUS                     → 16c /  96GB / 48h
```

---

## Routing rules (summary)

| Assay | Modality | Route |
|---|---|---|
| CITE, GEX | GEX, ADT, HTO, VDJ-*, CRISPR | CELLRANGER_MULTI → QC_GEX |
| DOGMA | GEX, ADT, HTO | CELLRANGER_MULTI → QC_GEX |
| DOGMA | ATAC | CELLRANGER_ATAC → QC_ATAC |
| Multiome | GEX, ADT | CELLRANGER_MULTI → QC_GEX |
| Multiome | ATAC | CELLRANGER_ATAC → QC_ATAC |
| ATAC | ATAC | CELLRANGER_ATAC → QC_ATAC |
| ASAP | ATAC | CELLRANGER_ATAC → QC_ATAC + COUNT_ADT |
| ASAP | ADT, HTO | COUNT_ADT (kallisto) |
| Any | GENO | skipped |

**DOGMA ADT never goes to kallisto.** It is handled by cellranger multi on the GEX run.

**One pipeline invocation = one BCL folder** (or equivalent). DOGMA GEX and DOGMA ATAC runs are separate invocations. Integration of GEX + ATAC outputs for DOGMA is done downstream in R.

---

## What NOT to do

- Do not add `run_id`, `bcl_path`, or any new column to the samplesheet — it must stay identical to `OSCAR/bash/metadata.csv` format
- Do not merge DOGMA GEX + ATAC outputs in Nextflow — independent runs, R integration downstream
- Do not route DOGMA ADT to kallisto — it goes to cellranger multi
- Do not invent Nextflow operators; confirmed operators: `.branch`, `.join`, `.groupTuple`, `.combine`, `.flatten`, `.multiMap`, `.mix`, `.collect`, `.map`, `.filter`, `.set`, `groupKey()`
- Every process must emit `versions.yml`
- `n_donors=NA` → treat as 1; do not run genotyping
- `adt_csv_path` is null when `--adt_files_dir` is not set; CELLRANGER_MULTI receives `file('NO_FILE')` sentinel and skips the `[feature]` section
- ATAC whitelist path is hardcoded in `nextflow.config`; do not try to extract it from the container at runtime
