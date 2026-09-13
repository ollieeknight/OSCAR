#!/usr/bin/env python3
"""Check the TotalSeq barcode CSVs the feature barcode generator reads.

BioLegend publishes no API, and its site blocks scripted requests, so these
files are refreshed by hand from the Barcode Look-Up tool:

    https://www.biolegend.com/en-us/totalseq/barcode-lookup

Pick a format tab, choose "download all items as an excel document", save it
as CSV with the columns below, then run this to check it before committing.

    python3 assets/check_totalseq_csvs.py             # check the committed files
    python3 assets/check_totalseq_csvs.py new_a.csv   # check a fresh download

Exits non-zero if anything is wrong. With --diff OLD NEW it reports what
changed between two versions of one format instead.
"""
import csv
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent / "docs-src/tools"
COLUMNS = ["catalogue_number", "totalseq_id", "marker", "clone", "reactivity",
           "barcode_sequence"]
BARCODE_LEN = 15


def read(path):
    # utf-8-sig: the BioLegend export carries a BOM.
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def check(path):
    """Return a list of problems. Empty means the file is good."""
    rows = read(path)
    bad = []
    if not rows:
        return [f"{path.name}: no rows"]

    # helpers/functions.js splits each line on ',' and unpacks by position, so
    # the column order is load-bearing and a comma inside a field shifts every
    # field after it.
    got = [c for c in rows[0].keys() if c]
    if got[:len(COLUMNS)] != COLUMNS:
        bad.append(f"{path.name}: columns are {got[:len(COLUMNS)]}, expected {COLUMNS}")

    seen_barcode, seen_cat = {}, {}
    for i, row in enumerate(rows, start=2):
        where = f"{path.name}:{i}"
        bc = (row.get("barcode_sequence") or "").strip()

        if len(bc) != BARCODE_LEN:
            bad.append(f"{where}: barcode {bc!r} is {len(bc)}bp, expected {BARCODE_LEN}")
        if set(bc) - set("ACGT"):
            bad.append(f"{where}: barcode {bc!r} has non-ACGT characters")

        # Two antibodies sharing a barcode cannot be told apart after counting.
        if bc in seen_barcode:
            bad.append(f"{where}: barcode {bc} already used on line {seen_barcode[bc]}")
        else:
            seen_barcode[bc] = i

        cat = (row.get("catalogue_number") or "").strip()
        if cat in seen_cat:
            bad.append(f"{where}: catalogue number {cat} already used on line {seen_cat[cat]}")
        else:
            seen_cat[cat] = i

        for col in COLUMNS:
            if not (row.get(col) or "").strip():
                bad.append(f"{where}: empty {col}")
            elif "," in (row.get(col) or ""):
                bad.append(f"{where}: {col} contains a comma, which breaks the generator")
    return bad


def diff(old_path, new_path):
    """Print what changed between two versions of one format's CSV."""
    old = {r["catalogue_number"]: r for r in read(old_path)}
    new = {r["catalogue_number"]: r for r in read(new_path)}

    added = sorted(set(new) - set(old))
    removed = sorted(set(old) - set(new))
    changed = [c for c in sorted(set(old) & set(new))
               if old[c]["barcode_sequence"] != new[c]["barcode_sequence"]]

    print(f"{len(old)} rows -> {len(new)} rows")
    for c in added:
        print(f"  added   {c}  {new[c]['marker']} ({new[c]['clone']})")
    for c in removed:
        print(f"  removed {c}  {old[c]['marker']} ({old[c]['clone']})")
    for c in changed:
        # A changed barcode silently miscounts every panel already using it.
        print(f"  BARCODE CHANGED {c} {new[c]['marker']}: "
              f"{old[c]['barcode_sequence']} -> {new[c]['barcode_sequence']}")
    if not (added or removed or changed):
        print("  no changes")


def main():
    args = sys.argv[1:]
    if args and args[0] == "--diff":
        if len(args) != 3:
            sys.exit("usage: check_totalseq_csvs.py --diff OLD.csv NEW.csv")
        diff(Path(args[1]), Path(args[2]))
        return 0

    paths = [Path(a) for a in args] or sorted(TOOLS.glob("totalseq_*.csv"))
    problems = []
    for p in paths:
        found = check(p)
        print(f"{p.name}: {len(read(p))} rows, {len(found)} problem(s)")
        problems += found

    for p in problems:
        print(f"  {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
