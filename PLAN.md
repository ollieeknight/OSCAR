# OSCAR Repo Restructure Plan

**Goal:** Promote `nextflow/` as the primary pipeline. Archive legacy bash/conda/apptainer/tests under `old/`.

---

## Current → Target Structure

```
BEFORE                          AFTER
──────────────────────────────  ─────────────────────────────────
OSCAR/                          OSCAR/
├── README.md                   ├── README.md           (updated)
├── apptainer/      ──────────► old/apptainer/
├── bash/           ──────────► old/bash/
├── conda/          ──────────► old/conda/
├── nextflow/                   ├── main.nf             (promoted)
│   ├── main.nf     ──────────► ├── nextflow.config     (promoted)
│   ├── nextflow.config         ├── CLAUDE.md           (promoted, paths updated)
│   ├── CLAUDE.md               ├── apptainer/          (from nextflow/apptainer/)
│   ├── apptainer/  ──────────► ├── assets/             (promoted)
│   ├── assets/     ──────────► ├── modules/            (promoted)
│   ├── modules/    ──────────► ├── subworkflows/       (promoted)
│   └── subworkflows/           ├── reference/          (unchanged)
├── reference/                  ├── site/               (unchanged)
├── site/                       ├── templates/          (unchanged)
├── templates/                  └── old/
└── tests/          ──────────►     ├── apptainer/
                                    ├── bash/
                                    ├── conda/
                                    └── tests/
```

---

## Phase 1 — Archive Legacy Folders

**What:** Move old bash-era folders into `old/` using `git mv` to preserve history.

**Steps:**
```bash
mkdir old
git mv apptainer  old/apptainer
git mv bash       old/bash
git mv conda      old/conda
git mv tests      old/tests
```

**Verification:**
- `ls old/` shows all four folders
- `ls` at root no longer shows `apptainer/`, `bash/`, `conda/`, `tests/`
- `git status` shows renames, no deletions

**Anti-patterns:**
- Do NOT use `mv` — use `git mv` to preserve history
- Do NOT delete; archive only

---

## Phase 2 — Promote nextflow/ to Root

**What:** Move all contents of `nextflow/` up to repo root using `git mv`.

**Steps (order matters — handle apptainer collision first):**
```bash
# 1. Promote nextflow/apptainer/ as the new root apptainer/
#    (old one is already in old/apptainer/ from Phase 1)
git mv nextflow/apptainer  apptainer

# 2. Promote remaining nextflow contents
git mv nextflow/assets       assets
git mv nextflow/modules      modules
git mv nextflow/subworkflows subworkflows
git mv nextflow/main.nf      main.nf
git mv nextflow/nextflow.config nextflow.config
git mv nextflow/CLAUDE.md    CLAUDE.md

# 3. Remove now-empty nextflow/ dir
#    (check it's empty first: ls nextflow/)
rmdir nextflow
```

**Files to NOT move** (generated/ephemeral, should be in .gitignore):
- `nextflow/.nextflow.log`
- `nextflow/.nextflow/` (cache directory if present)

**Verification:**
- `ls` at root shows `main.nf`, `nextflow.config`, `CLAUDE.md`, `apptainer/`, `assets/`, `modules/`, `subworkflows/`
- `ls nextflow/` fails (dir removed)
- `git status` shows renames only

---

## Phase 3 — Update CLAUDE.md

**What:** CLAUDE.md has a directory structure section that references `nextflow/` as a prefix. After promotion these paths are no longer correct.

**File:** `CLAUDE.md` (now at repo root)

**Section to update:** `## Directory structure` — the tree currently shows:
```
nextflow/
├── main.nf
├── nextflow.config
...
```
Change to reflect actual root layout (no `nextflow/` prefix).

**Also check and update:**
- Any relative path references like `../bash/` in the "Source bash scripts live in..." line — update to note bash scripts are archived in `old/bash/`
- The apptainer subdirectory path comment: `apptainer/mgatk2_docker/` stays correct

**Verification:**
- All paths in CLAUDE.md directory tree match actual file locations
- `grep "nextflow/" CLAUDE.md` returns zero results (no stale path prefixes)

---

## Phase 4 — Update README.md

**What:** Root `README.md` likely references bash scripts as primary usage. Update to reflect Nextflow as primary.

**File:** `README.md` (repo root)

**Changes needed:**
- Lead with Nextflow invocation as primary usage
- Add note that bash scripts are archived in `old/` for reference only
- Update any directory structure diagram if present

**Verification:**
- README accurately describes current repo structure
- No references to `bash/`, `conda/`, `apptainer/` at root level

---

## Phase 5 — Check .gitignore

**What:** Ensure generated Nextflow files are gitignored.

**File:** `.gitignore` at repo root (create if absent)

**Entries to add if missing:**
```
.nextflow/
.nextflow.log
work/
results/
```

**Verification:**
- `git status` does not show `.nextflow.log` or `work/` as untracked

---

## Phase 6 — Final Verification

```bash
# Structure check
ls                          # main.nf, nextflow.config, CLAUDE.md, modules/, subworkflows/, assets/, apptainer/, old/, reference/, site/, templates/
ls old/                     # apptainer/, bash/, conda/, tests/

# No stale references in CLAUDE.md
grep -r "nextflow/modules"  CLAUDE.md    # should be empty
grep -r "nextflow/main"     CLAUDE.md    # should be empty

# git status clean (only renames + edits, no unexpected deletions)
git status
git diff --stat HEAD
```

fix these warnings, we will never used these containers with the nextflow pipeline:
WARN: Access to undefined parameter `container_asap` -- Initialise it to a default value eg. `params.container_asap = some_value`
WARN: Access to undefined parameter `container_cellbender` -- Initialise it to a default value eg. `params.container_cellbender = some_value`

---

## Notes

- `reference/`, `site/`, `templates/` are untouched — they remain at root and are still relevant
- The `old/` folder name is intentional — clear that it's archived, not deleted
- Nextflow invocation command `nextflow run main.nf` stays identical after promotion (was `nextflow run nextflow/main.nf` before)
