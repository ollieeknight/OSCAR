#! /bin/bash

# Check if conda is installed and available in the PATH
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda and add it to your PATH."
    exit 1
fi

# Exit immediately if a command exits with a non-zero status
set -e

# Create a new conda environment named 'oscar_count' with specified packages from the 'bih-cubi' channel
conda create -y -n oscar_count -c https://repo.prefix.dev/romitools \
    bcl2fastq2 \
    fastqc \
    multiqc \
    bustools \
    pandas

# Activate the 'oscar_count' environment
conda activate oscar_count

# Install the 'bio' package using pip
pip install bio

# Export the 'oscar_count' environment to a YAML file
conda env export > oscar_count.yml

conda deactivate

# Create a new conda environment named 'oscar_qc' with specified packages from the 'nvidia' channel
conda create -y -n oscar_qc -c nvidia \
    python=3.7 \
    cellsnp-lite \
    numpy=1.19 \
    pandas \
    scipy \
    statsmodels \
    openjdk \
    r-base \
    r-data.table \
    r-matrix \
    bioconductor-genomicranges \
    bioconductor-summarizedexperiment \
    cuda-toolkit \
    cuda-nvcc

# Activate the 'oscar_qc' environment
conda activate oscar_qc

# Install additional packages using pip
pip install vireoSNP

# Export the 'oscar_qc' environment to a YAML file
conda env export > oscar_qc.yml

# Deactivate the current conda environment
conda deactivate