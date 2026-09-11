#!/usr/bin/env python3
"""Self-checks for bin/demux_qc.py. No cluster, no containers, no BCL data."""

import csv
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import demux_qc  # noqa: E402

REPO = os.path.join(os.path.dirname(__file__), "..")


def _tmp(name, text):
    path = os.path.join(TMP, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


TMP = tempfile.mkdtemp()

# ── parsing ───────────────────────────────────────────────────────────────
stats = _tmp(
    "Demultiplex_Stats.csv",
    "Lane,SampleID,Index,# Reads,# Perfect Index Reads\n"
    "1,GEX_A,GTAACATGCG-AGTGTTACCT,1000000,999000\n"
    "1,GEX_B,GTGGATCGCA-ACACTAAGGC,900000,899000\n"
    "1,ADT_A,ATTCAGAA,800000,799000\n"
    "1,HTO_A,TTCTGAAT,5000,4900\n"
    "1,Undetermined,,300000,0\n"
    "2,GEX_A,GTAACATGCG-AGTGTTACCT,500000,499000\n"
    "2,Undetermined,,10000,0\n",
)

per_lane = demux_qc.read_demultiplex_stats(stats)
assert set(per_lane) == {"1", "2"}, per_lane
assert per_lane["1"]["GEX_A"] == 1000000
assert per_lane["1"]["Undetermined"] == 300000
assert per_lane["2"]["GEX_A"] == 500000
print("OK: Demultiplex_Stats.csv parsed per lane")

# Percentages are taken against the lane total *including* Undetermined.
summary = demux_qc.summarise(per_lane)
lane1 = [r for r in summary if r["lane"] == "1"]
# 4-decimal rounding, so allow a hair of slack rather than demanding exactly 100.
assert abs(sum(float(r["percent_of_lane"]) for r in lane1) - 100.0) < 0.001
assert [r["library"] for r in lane1][0] == "GEX_A", "rows must be sorted by reads desc"
undet = [r for r in lane1 if r["library"] == "Undetermined"][0]
assert abs(float(undet["percent_of_lane"]) - 100.0 * 300000 / 3005000) < 1e-4
print("OK: summary percentages include Undetermined and sum to 100")

# A sample split across several index rows accumulates instead of overwriting.
multi = _tmp(
    "multi.csv",
    "Lane,SampleID,Index,# Reads\n1,SI_GA_A,AAACGGCG,10\n1,SI_GA_A,CCTACCAT,20\n",
)
assert demux_qc.read_demultiplex_stats(multi)["1"]["SI_GA_A"] == 30
print("OK: multi-row (4-sequence single-index) samples accumulate")

# ── dropout detection ─────────────────────────────────────────────────────
drops = demux_qc.find_dropouts(per_lane)
assert len(drops) == 1, drops
assert drops[0]["library"] == "HTO_A" and drops[0]["lane"] == "1"
print("OK: dropout flagged >100x below lane median")

# Undetermined must not be treated as a library, nor skew the median.
assert all(d["library"] != "Undetermined" for d in drops)

# Zero-read library is flagged without dividing by zero.
zero = demux_qc.read_demultiplex_stats(
    _tmp(
        "zero.csv",
        "Lane,SampleID,# Reads\n1,A,1000000\n1,B,900000\n1,C,0\n1,Undetermined,5\n",
    )
)
zdrops = demux_qc.find_dropouts(zero)
assert [d["library"] for d in zdrops] == ["C"], zdrops
assert "0 reads" in zdrops[0]["detail"]
print("OK: zero-read library flagged without ZeroDivisionError")

# A lane with a single declared library has no siblings -> no conclusion drawn.
solo = demux_qc.read_demultiplex_stats(
    _tmp("solo.csv", "Lane,SampleID,# Reads\n1,A,3\n1,Undetermined,9000000\n")
)
assert demux_qc.find_dropouts(solo) == []
print("OK: single-library lane produces no dropout warning")

# Evenly loaded lane is clean.
even = demux_qc.read_demultiplex_stats(
    _tmp("even.csv", "Lane,SampleID,# Reads\n1,A,100\n1,B,90\n1,C,80\n")
)
assert demux_qc.find_dropouts(even) == []
print("OK: evenly loaded lane produces no warnings")

# ── kit index lookup ──────────────────────────────────────────────────────
kits = demux_qc.load_kit_indexes(os.path.join(REPO, "assets", "indexes"))
assert kits.get("GTAACATGCG") == "SI-TT-A1", kits.get("GTAACATGCG")
assert kits.get("ATTACTCG") == "D701", kits.get("ATTACTCG")
assert "INDEX_NAME" not in kits and "INDEX(I7)" not in kits, "header row leaked"
print(f"OK: {len(kits)} kit index sequences loaded from assets/indexes")

# ── undeclared index detection ────────────────────────────────────────────
unknown = _tmp(
    "Top_Unknown_Barcodes.csv",
    "Lane,index,index2,# Reads\n"
    "1,GAATTCGT,,200000\n"      # real kit index, 6.6% of lane -> flag
    "1,ACGTACGT,,150000\n"      # 5% of lane but not a kit index -> ignore
    "1,TAATGCGC,,500\n",        # kit index but 0.017% of lane -> below threshold
)
und = demux_qc.find_undeclared(
    demux_qc.read_top_unknown(unknown), per_lane, kits, demux_qc.DEFAULT_UNKNOWN_PCT
)
assert [w["index"] for w in und] == ["GAATTCGT"], und
assert und[0]["reads"] == 200000
print("OK: undeclared kit index flagged; non-kit and sub-threshold ignored")

# ── concatenation of several groups' Reports ──────────────────────────────
# Each demux group writes its own Demultiplex_Stats.csv with the same filename,
# so DEMUX_QC stages them into numbered subdirs and concatenates. Only the
# first header survives; a leaked header row would parse as a bogus sample and
# a whole missing group would look like a normal small pool if not merged.
g1 = "Lane,SampleID,# Reads\n1,GEX_A,1000000\n1,Undetermined,5000\n"
g2 = "Lane,SampleID,# Reads\n1,ADT_A,900000\n1,Undetermined,4000\n"
merged = g1 + "".join(g2.splitlines(keepends=True)[1:])
per_lane_merged = demux_qc.read_demultiplex_stats(_tmp("merged.csv", merged))
assert set(per_lane_merged["1"]) == {"GEX_A", "ADT_A", "Undetermined"}, per_lane_merged
# Undetermined appears once per group and must accumulate, not overwrite.
assert per_lane_merged["1"]["Undetermined"] == 9000
assert "SampleID" not in per_lane_merged["1"], "duplicate header row leaked in as a sample"
print("OK: multi-group Reports concatenate, Undetermined accumulates, no header leak")

print("\nAll demux_qc self-checks passed.")
