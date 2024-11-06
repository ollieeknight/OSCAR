#! /bin/bash

# Check if conda is installed and available in the PATH
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda and add it to your PATH."
    exit 1
fi

# Exit immediately if a command exits with a non-zero status
set -e

# Create a new conda environment named 'oscar_count' with specified packages from the 'bih-cubi' channel
conda create -y -n oscar_count -c bih-cubi \
    bcl2fastq2 \
    fastqc \
    multiqc \
    bustools

# Activate the 'oscar_count' environment
conda activate oscar_count

# Install the 'bio' package using pip
pip install bio

# Create a directory named 'conda' if it doesn't exist
mkdir -p conda

# Export the 'oscar_count' environment to a YAML file
conda env export > conda/oscar_count.yml

# Create a new conda environment named 'oscar_qc' with specified packages from the 'nvidia' channel
conda create -y -n oscar_qc -c nvidia \
    python=3.7 \
    cellsnp-lite \
    numpy=1.19 \
    pandas \
    scipy \
    statsmodels \
    openjdk \
    r-base=4.2.3 \
    r-data.table \
    r-matrix \
    bioconductor-genomicranges \
    bioconductor-summarizedexperiment \
    cuda-toolkit \
    cuda-nvcc

# Activate the 'oscar_qc' environment
conda activate oscar_qc

# Install additional packages using pip
pip install vireoSNP cellbender mgatk

# Create a directory named 'conda' if it doesn't exist (already created earlier, but this ensures it exists)
mkdir -p conda

# Export the 'oscar_qc' environment to a YAML file
conda env export > conda/oscar_qc.yml

# Deactivate the current conda environment
conda deactivate
