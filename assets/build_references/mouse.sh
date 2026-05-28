#!/bin/bash
set -euo pipefail

BIND_PATHS="/charite-store-f/f-cc12-ag-romagnani,/sc-projects/sc-proj-cc12-ag-romagnani,/sc-scratch/sc-scratch-cc12-ag-romagnani,/home/knighto"

BEDTOOLS_IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-bedtools_2.27.1--h077b44d_9.img"
SAMTOOLS_IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-samtools-1.23.1--ha83d96e_0.img"
CELLRANGER_ARC_IMG="/sc-scratch/sc-scratch-cc12-ag-romagnani/apptainer_cache/quay.io-biocontainers-cellranger-arc_2.0.2.img"

OUTDIR="/sc-projects/sc-proj-cc12-ag-romagnani/ref/mm"
GENOME="GRCm38-hardmasked-optimised-arc"
BUILD="${OUTDIR}/${GENOME}-build"
SOURCE="${OUTDIR}/${GENOME}-source"
FASTA_NAME="Mus_musculus.GRCm38.dna_sm.primary_assembly"
FASTA_URL="https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/${FASTA_NAME}.fa.gz"
GTF_URL="https://storage.googleapis.com/generecovery/mouse_mm10_optimized_annotation_v2.gtf.gz"
BLACKLIST_URL="https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/mm10.full.blacklist.bed"
# MOTIFS_URL="https://testjaspar.uio.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

mkdir -p "${OUTDIR}" "${BUILD}" "${SOURCE}"

FASTA_IN="${SOURCE}/${FASTA_NAME}.fa"
if [ ! -f "${FASTA_IN}" ]; then
    echo "Downloading FASTA..."
    curl -sS "${FASTA_URL}" | zcat > "${FASTA_IN}"
fi

GTF_IN="${SOURCE}/mouse_mm10_optimized_annotation_v2.gtf"
if [ ! -f "${GTF_IN}" ]; then
    echo "Downloading GTF..."
    curl -sS "${GTF_URL}" | zcat > "${GTF_IN}"
fi

BLACKLIST_IN="${SOURCE}/mm10.full.blacklist.bed"
if [ ! -f "${BLACKLIST_IN}" ]; then
    echo "Downloading blacklist BED file..."
    curl -sS "${BLACKLIST_URL}" > "${BLACKLIST_IN}"
fi

FASTA_MOD="${BUILD}/${FASTA_NAME}.fa.mod"
if [ ! -f "${FASTA_MOD}" ]; then
    echo "Reformatting FASTA headers..."
    sed -E \
        's/^>(\S+).*/>\1 \1/; s/^>([0-9]+|[XY]) />chr\1 /; s/^>MT />chrM /' \
        "${FASTA_IN}" > "${FASTA_MOD}"
fi

FASTA_MASKED="${BUILD}/${FASTA_NAME}_hardmasked.fa.mod"
if [ ! -f "${FASTA_MASKED}" ]; then
    apptainer exec -B ${BIND_PATHS} "${BEDTOOLS_IMG}" \
        bedtools maskfasta \
            -fi "${FASTA_MOD}" \
            -bed "${BLACKLIST_IN}" \
            -fo "${FASTA_MASKED}"
fi

CONFIG_IN="${BUILD}/genome.config"
cat > "${CONFIG_IN}" <<CONFIGEOF
{
    organism: "Mus_musculus"
    genome: ["${GENOME}"]
    input_fasta: ["${FASTA_MASKED}"]
    input_gtf: ["${GTF_IN}"]
    # input_motifs: ""
    non_nuclear_contigs: ["chrM"]
}
CONFIGEOF

cd "${OUTDIR}"
apptainer exec -B ${BIND_PATHS} "${CELLRANGER_ARC_IMG}" \
    cellranger-arc mkref \
        --ref-version 'A' \
        --config "${CONFIG_IN}"

rm -r "${SOURCE}" "${BUILD}"
echo "Done. Reference built at ${OUTDIR}/${GENOME}/"