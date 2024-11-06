#!/bin/bash

# Check if the required tar.gz files are present in the parent directory
if [[ ! -f ../cellranger-8.0.0.tar.gz ]] && [[ ! -f ../cellranger-atac-2.1.0.tar.gz ]]; then
    # If not, prompt the user to download them
    echo "First, you need to download cellranger-8.0.0.tar.gz, cellranger-atac-2.1.0.tar.gz into this folder, then come back."
    echo "You can download them from here:"
    echo "https://www.10xgenomics.com/support/software/cell-ranger/downloads"
    echo "https://support.10xgenomics.com/single-cell-atac/software/downloads/latest"
fi

# Ask the user if they want to build the counting file
echo "Do you want to build the counting file? (Y/N)"
read -r choice

# Validate the user's input
while [[ ! $choice =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter Y or N."
    read -r choice
done

# Process the user's choice for building the counting file
if [[ "$choice" = "Y" ]] || [[ "$choice" = "y" ]]; then
    # If yes, build the counting file using apptainer
    echo "apptainer build oscar-count.sif recipe_oscar_count.sif"
    apptainer build oscar-count.sif recipe_oscar_count.sif
    echo ""
elif [[ "$choice" = "N" ]] || [[ "$choice" = "n" ]]; then
    # If no, do nothing
    :
fi

# Ask the user if they want to build the QC file
echo "Do you want to build the QC file? (Y/N)"
read -r choice

# Validate the user's input
while [[ ! $choice =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter Y or N."
    read -r choice
done

# Process the user's choice for building the QC file
if [[ "$choice" = "Y" ]] || [[ "$choice" = "y" ]]; then
    # If yes, build the QC file using apptainer
    echo ""
    echo "apptainer build oscar-qc.sif recipe_oscar_qc.sif"
    apptainer build oscar-qc.sif recipe_oscar_qc.sif
elif [[ "$choice" = "N" ]] || [[ "$choice" = "n" ]]; then
    # If no, do nothing
    :
fi
