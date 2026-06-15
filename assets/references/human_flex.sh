#!/bin/bash
# Build a cellranger (not cellranger-arc) reference for Flex runs that includes
# human GRCh38 + viral genomes (EBV/CMV/HHV6B/HHV7) needed for custom probe set.
#
# Follows the same pattern as build_references/human.sh but uses cellranger mkref
# and appends viral sequences before building.
#
# Usage:
#   bash build_flex_reference_with_viral.sh
#
# Output: ${OUTDIR}/GRCh38-hardmasked-flex-viral/
# Approx runtime: 2-3 h on BIH cluster

set -euo pipefail

BIND_PATHS="/charite-store-f/f-cc12-ag-romagnani,/sc-projects/sc-proj-cc12-ag-romagnani,/sc-scratch/sc-scratch-cc12-ag-romagnani,/home/knighto"

BEDTOOLS_IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-bedtools_2.27.1--h077b44d_9.img"
CELLRANGER_IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-nf-core-cellranger-10.0.0.img"

OUTDIR="/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs"
GENOME="GRCh38-hardmasked-flex-viral"
BUILD="${OUTDIR}/${GENOME}-build"
SOURCE="${OUTDIR}/${GENOME}-source"

FASTA_NAME="Homo_sapiens.GRCh38.dna_sm.primary_assembly"
FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${FASTA_NAME}.fa.gz"
GTF_URL="https://storage.googleapis.com/generecovery/human_GRCh38_optimized_annotation_v2.gtf.gz"
BLACKLIST_URL="https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed"

mkdir -p "${OUTDIR}" "${BUILD}" "${SOURCE}"

# ── Human FASTA ───────────────────────────────────────────────────────────────

FASTA_IN="${SOURCE}/${FASTA_NAME}.fa"
if [ ! -f "${FASTA_IN}" ]; then
    echo "[1/8] Downloading human FASTA..."
    curl -sS "${FASTA_URL}" | zcat > "${FASTA_IN}"
fi

GTF_IN="${SOURCE}/human_GRCh38_optimized_annotation_v2.gtf"
if [ ! -f "${GTF_IN}" ]; then
    echo "[2/8] Downloading human GTF..."
    curl -sS "${GTF_URL}" | zcat > "${GTF_IN}"
fi

BLACKLIST_IN="${SOURCE}/hg38.full.blacklist.bed"
if [ ! -f "${BLACKLIST_IN}" ]; then
    echo "[3/8] Downloading blacklist BED..."
    curl -sS "${BLACKLIST_URL}" > "${BLACKLIST_IN}"
fi

FASTA_MOD="${BUILD}/${FASTA_NAME}.fa.mod"
if [ ! -f "${FASTA_MOD}" ]; then
    echo "[4/8] Reformatting FASTA headers..."
    sed -E \
        's/^>(\S+).*/>\1 \1/; s/^>([0-9]+|[XY]) />chr\1 /; s/^>MT />chrM /' \
        "${FASTA_IN}" > "${FASTA_MOD}"
fi

FASTA_MASKED="${BUILD}/${FASTA_NAME}_hardmasked.fa.mod"
if [ ! -f "${FASTA_MASKED}" ]; then
    echo "[5/8] Masking blacklist regions..."
    apptainer exec -B "${BIND_PATHS}" "${BEDTOOLS_IMG}" \
        bedtools maskfasta \
            -fi "${FASTA_MOD}" \
            -bed "${BLACKLIST_IN}" \
            -fo "${FASTA_MASKED}"
fi

# ── Viral genomes ─────────────────────────────────────────────────────────────
# Accessions:
#   EBV   (HHV-4) : NC_007605.2   Epstein-Barr virus type 1
#   CMV   (HHV-5) : NC_006273.4   Human cytomegalovirus strain Merlin
#   HHV6B (HHV-6B): NC_000898.1   Human herpesvirus 6B strain Z29
#   HHV7  (HHV-7) : NC_001716.3   Human herpesvirus 7 strain JI
#
# Pre-download recommended (NCBI rate-limits unauthenticated requests):
#   mkdir -p ${VIRAL_PREFETCH_DIR}
#cd ${VIRAL_PREFETCH_DIR}
#for ACC in NC_007605.2 NC_006273.4 NC_000898.1 NC_001716.3; do
#curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" -o ${ACC}.fa
#sleep 3
#done
#mv NC_007605.2.fa chrEBV.fa
#mv NC_006273.4.fa chrCMV.fa
#mv NC_000898.1.fa chrHHV6B.fa
#mv NC_001716.3.fa chrHHV7.fa
#
# Script checks VIRAL_PREFETCH_DIR first; falls back to live NCBI download if absent.

VIRAL_PREFETCH_DIR="/sc-projects/sc-proj-cc12-ag-romagnani/ref/hs/viral_source"

VIRAL_FASTA="${BUILD}/viral_genomes.fa"
# Validate existing file: must have exactly 4 viral chromosome headers.
# A stale empty file from a previously failed run would pass -f but fail this check.
_viral_ok=false
if [ -f "${VIRAL_FASTA}" ] && [ "$(grep -c '^>' "${VIRAL_FASTA}" 2>/dev/null || echo 0)" -eq 4 ]; then
    _viral_ok=true
fi
if [ "${_viral_ok}" = "false" ]; then
    [ -f "${VIRAL_FASTA}" ] && echo "  WARNING: existing ${VIRAL_FASTA} is invalid (removing and regenerating)"
    rm -f "${VIRAL_FASTA}"
    echo "[6/8] Building viral genome FASTA..."

    declare -A VIRAL_ACC=(
        ["chrEBV"]="NC_007605.2"
        ["chrCMV"]="NC_006273.4"
        ["chrHHV6B"]="NC_000898.1"
        ["chrHHV7"]="NC_001716.3"
    )

    > "${VIRAL_FASTA}"
    for CHR in chrEBV chrCMV chrHHV6B chrHHV7; do
        ACC="${VIRAL_ACC[$CHR]}"
        # Accept either chrXXX.fa or NC_accession.fa naming
        PRE_CHR="${VIRAL_PREFETCH_DIR}/${CHR}.fa"
        PRE_ACC="${VIRAL_PREFETCH_DIR}/${ACC}.fa"

        if [ -f "${PRE_CHR}" ]; then
            echo "  Using pre-downloaded ${PRE_CHR}"
            SRC="${PRE_CHR}"
        elif [ -f "${PRE_ACC}" ]; then
            echo "  Using pre-downloaded ${PRE_ACC}"
            SRC="${PRE_ACC}"
        else
            echo "  Fetching ${ACC} from NCBI → ${CHR}"
            TMP_FA="${BUILD}/tmp_${CHR}.fa"
            curl --fail --retry 3 --retry-delay 10 -sS \
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" \
                -o "${TMP_FA}"
            SRC="${TMP_FA}"
            sleep 2
        fi

        if ! head -1 "${SRC}" | grep -q "^>"; then
            echo "ERROR: ${SRC} is not valid FASTA (first line: $(head -1 "${SRC}"))"
            exit 1
        fi

        # Re-header to chrXXX and strip blank lines
        awk -v chr="${CHR}" '/^>/{print ">" chr " " chr; next} /^$/{next} {print}' \
            "${SRC}" >> "${VIRAL_FASTA}"

        rm -f "${BUILD}/tmp_${CHR}.fa"
    done

    N_HEADERS=$(grep -c "^>" "${VIRAL_FASTA}")
    if [ "${N_HEADERS}" -ne 4 ]; then
        echo "ERROR: Expected 4 viral chromosomes, got ${N_HEADERS}:"
        grep "^>" "${VIRAL_FASTA}"
        exit 1
    fi
    echo "  Viral FASTA OK: ${N_HEADERS} chromosomes"
fi

# ── Viral GTF ─────────────────────────────────────────────────────────────────
# Minimal gene models — one gene / one transcript / one exon spanning the entire
# chromosome. Coordinates are fake (1 to chrom_end) but cellranger only needs
# the gene_id to exist for Flex probe counting; exact coordinates do not matter
# because probe reads are assigned by barcode, not by alignment position.
#
# gene_id values MUST match gene_id column in custom_probe_example.csv exactly.
# IMPORTANT: if you update the probe set, update this GTF section too.

VIRAL_GTF="${BUILD}/viral_genes.gtf"
if [ ! -f "${VIRAL_GTF}" ]; then
    echo "[7/8] Generating viral gene GTF entries..."

    # Get chromosome lengths from the viral FASTA we just downloaded.
    # Pass paths as positional args so the single-quoted heredoc stays unexpanded.
    python3 - "${VIRAL_FASTA}" "${VIRAL_GTF}" << 'PYEOF'
import sys

viral_fasta, viral_gtf = sys.argv[1], sys.argv[2]

chrom_lens = {}
current = None
length = 0
with open(viral_fasta) as fh:
    for line in fh:
        if line.startswith(">"):
            if current:
                chrom_lens[current] = length
            current = line.split()[0][1:]
            length = 0
        else:
            length += len(line.strip())
if current:
    chrom_lens[current] = length

# gene_id → (chromosome, gene_name)
# gene_id must exactly match probe CSV gene_id column
GENES = [
    # EBV (HHV-4) — chrEBV
    ("HHV4_BZLF1",  "chrEBV",   "HHV4_BZLF1"),
    ("HHV4_EBNA-2", "chrEBV",   "HHV4_EBNA-2"),
    ("HHV4_LMP-1",  "chrEBV",   "HHV4_LMP-1"),
    ("HHV4_BHLF1",  "chrEBV",   "HHV4_BHLF1"),
    # CMV (HHV-5) — chrCMV
    ("HHV5_gp108",  "chrCMV",   "HHV5_UL123"),
    ("HHV5_gp166",  "chrCMV",   "HHV5_RNA2.7"),
    ("HHV5_gp027",  "chrCMV",   "HHV5_UL22A"),
    ("HHV5_gp158",  "chrCMV",   "HHV5_US28"),
    ("HHV5_gp124",  "chrCMV",   "HHV5_UL138"),
    ("HHV5_gp041",  "chrCMV",   "HHV5_UL36"),
    # HHV-6B — chrHHV6B
    ("HHV6_gp091",  "chrHHV6B", "HHV6_U90"),
    ("HHV6_gp047",  "chrHHV6B", "HHV6_U38"),
    ("HHV6_gp096",  "chrHHV6B", "HHV6_U100"),
    # HHV-7 — chrHHV7
    ("HHV7_gp84",   "chrHHV7",  "HHV7_U100"),
    ("HHV7_gp86",   "chrHHV7",  "HHV7_DR6"),
    ("HHV7_gp09",   "chrHHV7",  "HHV7_U11"),
    ("HHV7_gp38",   "chrHHV7",  "HHV7_U38"),
]

with open(viral_gtf, "w") as f:
    for gene_id, chrom, gene_name in GENES:
        end = chrom_lens.get(chrom, 200000)
        tx_id = gene_id + "_T1"
        attrs_gene = f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_biotype "protein_coding";'
        attrs_tx   = f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
        f.write(f'{chrom}\tcustom\tgene\t1\t{end}\t.\t+\t.\t{attrs_gene}\n')
        f.write(f'{chrom}\tcustom\ttranscript\t1\t{end}\t.\t+\t.\t{attrs_tx}\n')
        f.write(f'{chrom}\tcustom\texon\t1\t{end}\t.\t+\t.\t{attrs_tx}\n')

print(f"Wrote {len(GENES)} viral gene models to {viral_gtf}")
PYEOF
fi

# ── Combine human + viral ─────────────────────────────────────────────────────

FASTA_COMBINED="${BUILD}/genome_combined.fa"
GTF_COMBINED="${BUILD}/genes_combined.gtf"

# Validate existing combined FASTA: must contain all 4 viral contigs.
# Guards against a stale combined file built before viral genomes were appended.
_combined_ok=false
if [ -f "${FASTA_COMBINED}" ]; then
    _missing=0
    for _chr in chrEBV chrCMV chrHHV6B chrHHV7; do
        grep -q "^>${_chr}" "${FASTA_COMBINED}" || _missing=$((_missing + 1))
    done
    [ "${_missing}" -eq 0 ] && _combined_ok=true
fi
if [ "${_combined_ok}" = "false" ]; then
    [ -f "${FASTA_COMBINED}" ] && echo "  WARNING: ${FASTA_COMBINED} missing viral contigs — regenerating"
    echo "[8/8] Concatenating human + viral FASTA and GTF..."
    cat "${FASTA_MASKED}" "${VIRAL_FASTA}" > "${FASTA_COMBINED}"
    cat "${GTF_IN}" "${VIRAL_GTF}" > "${GTF_COMBINED}"
fi

# ── cellranger mkref ──────────────────────────────────────────────────────────

cd "${OUTDIR}"
apptainer exec -B "${BIND_PATHS}" "${CELLRANGER_IMG}" \
    cellranger mkref \
        --genome="${GENOME}" \
        --fasta="${FASTA_COMBINED}" \
        --genes="${GTF_COMBINED}" \
        --ref-version="1.0" \
        --nthreads=16 \
        --memgb=64

rm -r "${SOURCE}" "${BUILD}"