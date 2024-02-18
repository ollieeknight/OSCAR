#!/bin/bash

if ! command -v conda &> /dev/null; then
    echo "Error: conda is not found in the PATH. Please ensure conda is installed and added to your PATH."
    exit 1
fi

conda create -y -n oscar_count -c bih-cubi bcl2fastq2 fastqc multiqc bustools
conda activate oscar_count
pip install bio
conda env export > conda/oscar_count.yml
conda deactivate

####

conda create -y -n oscar_qc -c nvidia python=3.7 cellsnp-lite numpy=1.19 pandas scipy statsmodels openjdk r-base=4.2.3 r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment cuda-toolkit cuda-nvcc
conda activate oscar_qc
pip install vireoSNP cellbender mgatk
conda env export > conda/oscar_qc.yml
conda deactivate
