#!/usr/bin/env python3
"""Generate docs-src/tools/helpers/schema.js from the pipeline's own definitions.

The samplesheet generator used to hard-code its dropdown options, so it drifted
from the pipeline and offered chemistries the pipeline rejects. Reading
lib/chemistry.nf and lib/samplesheet.nf keeps the two in step.

Run after changing valid assays, modalities, or chemistries:

    python3 assets/generate_tool_schema.py
"""

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "docs-src/tools/helpers/schema.js"

# Human-readable labels. A chemistry absent here falls back to its raw name,
# so a new chemistry still appears in the dropdown without touching this file.
CHEMISTRY_LABELS = {
    "SC3Pv2": "Single cell 3', v2",
    "SC3Pv3": "Single cell 3', v3",
    "SC3Pv4": "Single cell 3', v4 (GEM-X)",
    "SC5P": "Single cell 5', v2",
    "SC5P-R2": "Single cell 5', R2-only, v2",
    "SC5Pv3": "Single cell 5', v3",
    "ARCv1": "Multiome/DOGMA, v1",
    "ARC-v1": "Multiome/DOGMA, v1 (hyphenated)",
    "Flex-v2-R1": "Flex (Fixed RNA), v2, R1",
    "Flex-v2-RNA-R2": "Flex (Fixed RNA), v2, RNA-R2",
    "ATAC": "ATAC/ASAP: v1, v1.1 and v2",
    "NA": "Not applicable",
}

MODALITY_LABELS = {
    "GEX": "Gene expression",
    "ADT": "Surface antibody",
    "HTO": "Surface hashtag",
    "ATAC": "ATAC",
    "VDJ-T": "VDJ, T cell",
    "VDJ-B": "VDJ, B cell",
    "CRISPR": "CRISPR guide capture",
    "GENO": "GoTChA genotyping",
}


def groovy_list(path: Path, func: str) -> list[str]:
    """Pull ['a', 'b'] out of a one-line Groovy `def func() { [...] }`."""
    text = path.read_text()
    m = re.search(rf"def {func}\(\)\s*{{\s*(\[[^\]]*\])", text)
    if not m:
        sys.exit(f"ERROR: could not find {func}() in {path}")
    return re.findall(r"'([^']+)'", m.group(1))


def chemistry_keys(path: Path) -> list[str]:
    """Pull the top-level keys out of chemistry_registry()."""
    text = path.read_text()
    start = text.index("def chemistry_registry()")
    end = text.index("\n}", start)
    return re.findall(r"^\s{8}'([^']+)':\s*\[", text[start:end], re.MULTILINE)


def index_codes(assets: Path) -> list[str]:
    """Every 10x kit code shipped in assets/indexes/."""
    codes = []
    for csv in sorted(assets.glob("*_Index_Kit_*.csv")):
        for line in csv.read_text().splitlines():
            if not line.strip() or line.startswith("index_name"):
                continue
            codes.append(line.split(",")[0].strip())
    return codes


def main() -> None:
    lib = REPO / "lib"
    schema = {
        "assays": groovy_list(lib / "samplesheet.nf", "valid_assays"),
        "indexTypes": groovy_list(lib / "samplesheet.nf", "valid_index_types"),
        "columns": [],
        "modalities": [],
        "chemistries": [],
        "indexCodes": index_codes(REPO / "assets/indexes"),
    }

    for mod in groovy_list(lib / "samplesheet.nf", "valid_modalities"):
        schema["modalities"].append({"value": mod, "label": MODALITY_LABELS.get(mod, mod)})

    for chem in chemistry_keys(lib / "chemistry.nf"):
        schema["chemistries"].append(
            {"value": chem, "label": CHEMISTRY_LABELS.get(chem, chem)}
        )

    # Column order defines the CSV header, so it comes from the pipeline too.
    text = (lib / "samplesheet.nf").read_text()
    m = re.search(r"def samplesheet_required_columns\(\)\s*{\s*(\[.*?\])", text, re.S)
    if not m:
        sys.exit("ERROR: could not find samplesheet_required_columns()")
    schema["columns"] = re.findall(r"'([^']+)'", m.group(1))

    body = json.dumps(schema, indent=2)
    OUT.write_text(
        "// GENERATED FILE - do not edit.\n"
        "// Regenerate with: python3 assets/generate_tool_schema.py\n"
        "// Source: lib/samplesheet.nf, lib/chemistry.nf, assets/indexes/\n"
        f"const OSCAR_SCHEMA = {body};\n"
    )
    print(
        f"Wrote {OUT.relative_to(REPO)}: "
        f"{len(schema['assays'])} assays, "
        f"{len(schema['modalities'])} modalities, "
        f"{len(schema['chemistries'])} chemistries, "
        f"{len(schema['indexCodes'])} index codes"
    )


if __name__ == "__main__":
    main()
