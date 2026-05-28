# PLAN: Fix SCRUBLET Process — Replace h5py with scanpy

## Problem

`QC_GEX:SCRUBLET` fails with `ModuleNotFoundError: No module named 'h5py'`.
The `quay.io/biocontainers/scrublet:0.2.3` container lacks h5py.
The process reads CellBender `output.h5` via h5py, then runs the scrublet Python API.

## Solution

Switch to the scanpy community container (which bundles h5py + scrublet) and rewrite
the inline Python to use `scanpy.read_10x_h5()` + `sc.pp.scrublet()`.
No interface changes: input is still the CellBender h5; output is still `doublets.csv`
with columns `barcode` (index), `doublet_score`, `is_gex_doublet`.

**Reference implementation**: nf-core/scdownstream `SCANPY_SCRUBLET` module (provided by user).
- Container: `community.wave.seqera.io/library/python_pyyaml_scanpy_scikit-image:750e7b74b6d036e4`
- Pattern: set `MPLCONFIGDIR` + `NUMBA_CACHE_DIR` before scanpy import (apptainer requirement)
- `sc.pp.scrublet()` adds `doublet_score` and `predicted_doublet` to `adata.obs`

---

## Phase 0: Documentation Discovery (COMPLETE)

**Findings confirmed from nf-core module and user-provided error:**

| Item | Detail |
|------|--------|
| Broken container | `quay.io/biocontainers/scrublet:0.2.3--py38h0213d0e_0` — no h5py |
| Replacement container | `community.wave.seqera.io/library/python_pyyaml_scanpy_scikit-image:750e7b74b6d036e4` |
| Read CellBender h5 | `sc.read_10x_h5(filename)` — CellBender uses 10x sparse format under `matrix/` group |
| Run scrublet | `sc.pp.scrublet(adata, expected_doublet_rate=0.08)` |
| Doublet score col | `adata.obs['doublet_score']` (float) |
| Doublet call col | `adata.obs['predicted_doublet']` (bool) → rename to `is_gex_doublet` |
| Apptainer fix | Export `MPLCONFIGDIR=./tmp/mpl` and `NUMBA_CACHE_DIR=./tmp/numba` before import |
| Container param | `params.container_scrublet` in `nextflow.config:42` |

**Files to edit:**
- `modules/qc.nf` lines ~252–303 (SCRUBLET process)
- `nextflow.config` line 42 (container_scrublet value)

---

## Phase 1: Edit `modules/qc.nf` — SCRUBLET process

### What to change

Replace the entire script block of the SCRUBLET process with scanpy-based logic.

**Old container line** (referenced via param — no direct change needed here, param updated in Phase 2):
```
container "${params.container_scrublet}"
```

**Old script block** (h5py approach, broken):
```python
import h5py
with h5py.File('${cellbender_h5}', 'r') as f:
    ...
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.08)
doublet_scores, predicted_doublets = scrub.doublet_detector()
```

**New script block** (scanpy approach):
```groovy
    script:
    """
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    python << 'PYEOF'
import os
import scanpy as sc
import pandas as pd

adata = sc.read_10x_h5('${cellbender_h5}')
sc.pp.scrublet(adata, expected_doublet_rate=0.08)

df = pd.DataFrame({
    'doublet_score': adata.obs['doublet_score'],
    'is_gex_doublet': adata.obs['predicted_doublet'].astype(bool)
}, index=adata.obs_names)
df.index.name = 'barcode'
df.to_csv('doublets.csv')
PYEOF

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)" 2>&1)
END_VERSIONS
    """
```

### Verification checklist

- [ ] `modules/qc.nf` SCRUBLET script block no longer imports `h5py` or `scrublet`
- [ ] Script exports `MPLCONFIGDIR` and `NUMBA_CACHE_DIR` before python invocation
- [ ] Output CSV declaration unchanged: `path("doublets.csv"), emit: doublets`
- [ ] versions.yml reports scanpy version (not scrublet)

---

## Phase 2: Update `nextflow.config` — container param

### What to change

Line 42:
```
# Old
container_scrublet     = "quay.io/biocontainers/scrublet:0.2.3--py38h0213d0e_0"

# New
container_scrublet     = "community.wave.seqera.io/library/python_pyyaml_scanpy_scikit-image:750e7b74b6d036e4"
```

The param name `container_scrublet` stays the same — no downstream references to update.

### Verification checklist

- [ ] `nextflow.config` line 42 points to new community container
- [ ] `grep container_scrublet nextflow.config` shows new URL

---

## Phase 3: Verify end-to-end

### Test approach

Re-run a single SCRUBLET task via Nextflow `-resume` to pick up the failed task with new container.

```bash
# From pipeline launch dir, resume from last checkpoint
nextflow run main.nf -resume [original params] -entry OSCAR
```

Or manually test the Python in the work dir:

```bash
# In the failed work dir: /sc-scratch/.../nf_work_oscar/eb/7c21741.../
apptainer exec community.wave.seqera.io/... python << 'EOF'
import os; os.environ["MPLCONFIGDIR"]="./tmp/mpl"; os.environ["NUMBA_CACHE_DIR"]="./tmp/numba"
import scanpy as sc
adata = sc.read_10x_h5('output.h5')
sc.pp.scrublet(adata, expected_doublet_rate=0.08)
print(adata.obs[['doublet_score','predicted_doublet']].head())
EOF
```

### Verification checklist

- [ ] SCRUBLET process completes with exit status 0
- [ ] `doublets.csv` exists with columns: `barcode`, `doublet_score`, `is_gex_doublet`
- [ ] Row count matches cell count from CellBender output
- [ ] `QC_GEX:CELLSNP_LITE` and downstream processes proceed normally after resume

---

## Anti-patterns to avoid

- Do NOT use `scrub.doublet_detector()` — that's the old scrublet Python API, not available in scanpy
- Do NOT use `sc.read_h5ad()` — CellBender output is `.h5` (10x format), not `.h5ad`; use `sc.read_10x_h5()`
- Do NOT skip `MPLCONFIGDIR`/`NUMBA_CACHE_DIR` exports — scanpy/numba will fail in apptainer without writable tmp dirs
- Do NOT rename `params.container_scrublet` — it would require updating `modules/qc.nf` container line too

---

## Notes

- The `doublets` channel is emitted in `subworkflows/qc_gex.nf:42` but not yet consumed in `main.nf`
  (output is published to disk only). No downstream channel consumers need updating.
- Container pull will happen automatically on first run via Nextflow's apptainer cache at
  `$NXF_APPTAINER_CACHEDIR` (already configured on the cluster).
