#!/bin/bash

# Define default paths and variables
CONDA_ENV_PATH="$HOME/work/bin/miniforge3/envs/genome_processing"
DEFAULT_CELLRANGER_PATH="$HOME/group/work/bin/cellranger-arc-2.0.2"
CELLRANGER_PATH=${1:-$DEFAULT_CELLRANGER_PATH}
GENOME="GRCh38-hardmasked-optimised-arc"
BUILD="${GENOME}-build"
SOURCE="${GENOME}-source"
FASTA_NAME="Homo_sapiens.GRCh38.dna_sm.primary_assembly"
FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${FASTA_NAME}.fa.gz"
GTF_URL="https://storage.googleapis.com/generecovery/human_GRCh38_optimized_annotation_v2.gtf.gz"
MOTIFS_URL="https://testjaspar.uio.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

# Create conda environment if it doesn't exist
if [ ! -d "$CONDA_ENV_PATH" ]; then
    conda create -y -n genome_processing bcftools samtools bedtools bwa
fi

# Check if conda environment was created
if [ ! -d "$CONDA_ENV_PATH" ]; then
    echo "Error: conda env 'genome_processing' still does not exist"
    exit 1
fi

# Activate the conda environment
conda activate genome_processing

# Check if Cell Ranger is installed
if [ ! -d "$CELLRANGER_PATH" ]; then
    echo "Please make sure Cell Ranger is properly pathed"
    exit 1
fi

# Update PATH for Cell Ranger
export PATH="$CELLRANGER_PATH:$PATH"

# Create directories for build and source
mkdir -p "${BUILD}" "${SOURCE}"

# Download FASTA, GTF, and motifs files if not already present
FASTA_IN="${SOURCE}/${FASTA_NAME}.fa"
GTF_IN="${SOURCE}/human_GRCh38_optimized_annotation_v2.gtf"
MOTIFS_IN="${SOURCE}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

if [ ! -f "${FASTA_IN}" ]; then
    curl -sS "${FASTA_URL}" | zcat > "${FASTA_IN}"
fi

if [ ! -f "${GTF_IN}" ]; then
    curl -sS "${GTF_URL}" | zcat > "${GTF_IN}"
fi

if [ ! -f "${MOTIFS_IN}" ]; then
    curl -sS "${MOTIFS_URL}" > "${MOTIFS_IN}"
fi

# Modify FASTA file
FASTA_MOD="${BUILD}/$(basename "${FASTA_IN}").mod"
sed -E 's/^>(\S+).*/>\1 \1/; s/^>([0-9]+|[XY]) />chr\1 /; s/^>MT />chrM /' "${FASTA_IN}" > "${FASTA_MOD}"

# Modify motifs file
MOTIFS_MOD="${BUILD}/$(basename "${MOTIFS_IN}").mod"
awk 'substr($1, 1, 1) == ">" { print ">" $2 "_" substr($1,2) } !substr($1, 1, 1) == ">" { print }' "${MOTIFS_IN}" > "${MOTIFS_MOD}"

# Download and mask the blacklist
curl -sS https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed > "${SOURCE}/hg38.full.blacklist.bed"
mv "${BUILD}/${FASTA_NAME}.fa.mod" "${BUILD}/${FASTA_NAME}_original.fa.mod"
bedtools maskfasta -fi "${BUILD}/${FASTA_NAME}_original.fa.mod" -bed "${SOURCE}/hg38.full.blacklist.bed" -fo "${BUILD}/${FASTA_NAME}_hardmasked.fa.mod"

# Create configuration file
CONFIG_IN="${BUILD}/genome.config"
cat <<EOF > "${CONFIG_IN}"
{
    organism: "Homo_sapiens"
    genome: ["${GENOME}"]
    input_fasta: ["${BUILD}/${FASTA_NAME}_hardmasked.fa.mod"]
    input_gtf: ["${GTF_IN}"]
    input_motifs: "${MOTIFS_MOD}"
    non_nuclear_contigs: ["chrM"]
}
EOF

# Run Cell Ranger to create reference
cellranger-arc mkref --ref-version 'A' --config "${CONFIG_IN}"

# Clean up
rm -r "${SOURCE}" "${BUILD}"
