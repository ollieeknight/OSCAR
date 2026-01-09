#!/bin/bash

# Set base directory and output directory
REFERENCE_FOLDER="/data/cephfs-1/home/users/knighto_c/group/work/ref/hs/"

mkdir -p ${REFERENCE_FOLDER}/source

cd ${REFERENCE_FOLDER}/source

wget "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
wget "https://storage.googleapis.com/generecovery/human_GRCh38_optimized_annotation_v2.gtf.gz"
wget "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
OUTPUT_DIR="${REFERENCE_FOLDER}/GRCh38_kallisto_optimised"

gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

export PATH="/data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/kallisto"

# Run kb ref command
kb ref --workflow=standard \
  -i "$OUTPUT_DIR/index.idx" \
  -g t2g.txt \
  -f1 Homo_sapiens.GRCh38.cdna.all.fa \
  --include-attribute gene_biotype:protein-coding \
  --include-attribute gene_biotype:lncRNA \
  --include-attribute gene_biotype:lincRNA \
  --include-attribute gene_biotype:antisense \
  --include-attribute gene_biotype:IG_LV_gene \
  --include-attribute gene_biotype:IG_V_gene \
  --include-attribute gene_biotype:IG_V_pseudogene \
  --include-attribute gene_biotype:IG_D_gene \
  --include-attribute gene_biotype:IG_J_gene \
  --include-attribute gene_biotype:IG_J_pseudogene \
  --include-attribute gene_biotype gene_biotype:IG_C_pseudogene \
  --include-attribute gene_biotype:TR_V_gene \
  --include-attribute gene_biotype:TR_V_pseudogene \
  --include-attribute gene_biotype:TR_D_gene \
  --include-attribute gene_biotype:TR_J_gene \
  --include-attribute gene_biotype:TR_J_pseudogene \
  --include-attribute gene_biotype:TR_C_gene \
  Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  human_GRCh38_optimized_annotation_v2.gtf.gz

echo "kb ref indexing completed. Index files are in $OUTPUT_DIR"