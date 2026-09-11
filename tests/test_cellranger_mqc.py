#!/usr/bin/env python3
"""Self-checks for bin/cellranger_mqc.py. No cluster, no cellranger required."""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import cellranger_mqc  # noqa: E402

TMP = tempfile.mkdtemp()

# Real `cellranger multi` layout and long-format CSV:
#   <library>/outs/per_sample_outs/<sample>/metrics_summary.csv
# Columns: Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value
HEADER = "Category,Library Type,Grouped By,Group Name,Metric Name,Metric Value\n"


def make_sample(library, sample, cells, mean_reads, extra=""):
    d = os.path.join(TMP, library, "outs", "per_sample_outs", sample)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "metrics_summary.csv"), "w", encoding="utf-8") as fh:
        fh.write(HEADER)
        fh.write(f'Cells,Gene Expression,Physical library ID,GEX,Cells,"{cells}"\n')
        fh.write(
            f'Cells,Gene Expression,Physical library ID,GEX,Mean reads per cell,"{mean_reads}"\n'
        )
        fh.write("Cells,Gene Expression,Physical library ID,GEX,Median genes per cell,1750\n")
        # Library-level rows: same metric namespace, different Category.
        fh.write('Library,Gene Expression,Physical library ID,GEX,Number of reads,"300,000,000"\n')
        fh.write("Library,Gene Expression,Physical library ID,GEX,Valid barcodes,97.3%\n")
        # Other library types must not leak into the Gene Expression view.
        fh.write("Cells,Antibody Capture,Physical library ID,ADT,Cells,99999\n")
        fh.write("Library,VDJ T,Physical library ID,VDJT,Number of reads,123\n")
        fh.write(extra)
    return d


make_sample("CITE_X_exp1_libA", "CITE_X_exp1_libA", "8,123", "45,000")
make_sample("CITE_X_exp1_libB", "CITE_X_exp1_libB", "7,004", "51,200")

# ── discovery ─────────────────────────────────────────────────────────────
entries = cellranger_mqc.find_metric_files(TMP)
assert len(entries) == 2, entries
assert sorted(s for s, _ in entries) == ["CITE_X_exp1_libA", "CITE_X_exp1_libB"]
print("OK: per_sample_outs/*/metrics_summary.csv discovered, sample id from dir name")

# ── parsing: Gene Expression only, Cells preferred over Library ───────────
m = cellranger_mqc.parse_metrics(
    os.path.join(TMP, "CITE_X_exp1_libA", "outs", "per_sample_outs",
                 "CITE_X_exp1_libA", "metrics_summary.csv")
)
assert m["Cells"] == "8,123", m["Cells"]
assert m["Cells"] != "99999", "Antibody Capture row leaked into Gene Expression metrics"
assert m["Median genes per cell"] == "1750"
# Library-level metric with no Cells-level counterpart is kept as a fallback.
assert m["Valid barcodes"] == "97.3%"
print("OK: Gene Expression rows only; Cells category wins, Library used as fallback")

# ── table build ───────────────────────────────────────────────────────────
present, rows = cellranger_mqc.build_table(entries)
assert set(rows) == {"CITE_X_exp1_libA", "CITE_X_exp1_libB"}
assert "Cells" in present and "Mean reads per cell" in present
# Only metrics actually present become columns.
assert "Sequencing saturation" not in present
# Column order follows WANTED, not file order.
assert present.index("Cells") < present.index("Median genes per cell")
print(f"OK: table built with {len(present)} metric column(s), ordered by WANTED")

# ── output file is valid MultiQC custom content ───────────────────────────
out = os.path.join(TMP, "cellranger_mqc.csv")
cellranger_mqc.write_mqc(out, present, rows)
text = open(out, encoding="utf-8").read()
assert text.startswith("# id: oscar_cellranger"), text[:40]
assert "# plot_type: 'table'" in text
assert "# section_name:" in text
body = [ln for ln in text.splitlines() if not ln.startswith("#")]
assert body[0].startswith("Sample,"), body[0]
assert len(body) == 3, body  # header + 2 samples
assert body[1].startswith("CITE_X_exp1_libA,"), body[1]
print("OK: MultiQC custom-content header and one row per sample")

# ── numeric coercion ──────────────────────────────────────────────────────
assert cellranger_mqc._num("1,234") == 1234.0
assert cellranger_mqc._num("92.5%") == 92.5
assert cellranger_mqc._num("N/A") is None
assert cellranger_mqc._num("") is None
print("OK: numeric coercion handles thousands separators, percents, N/A")

# ── empty tree is not an error ────────────────────────────────────────────
empty = tempfile.mkdtemp()
assert cellranger_mqc.find_metric_files(empty) == []
assert cellranger_mqc.main(["--root", empty, "--out", os.path.join(empty, "x.csv")]) == 0
assert not os.path.exists(os.path.join(empty, "x.csv")), "should not write an empty table"
print("OK: no cellranger output is a no-op, not a failure")

print("\nAll cellranger_mqc self-checks passed.")
