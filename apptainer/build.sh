#!/bin/bash/

if [[ ! -f ../cellranger-8.0.0.tar.gz && ! -f ../cellranger-atac-2.1.0.tar.gz ]]; then
    echo "First, you need to download cellranger-8.0.0.tar.gz, cellranger-atac-2.1.0.tar.gz into this folder, then come back."
    echo "You can download them from here:"
    echo "https://www.10xgenomics.com/support/software/cell-ranger/downloads"
    echo "https://support.10xgenomics.com/single-cell-atac/software/downloads/latest"
fi

echo "Do you want to build the counting file? (Y/N)"
read -r choice

while [[ ! $choice =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter Y or N."
    read -r choice
done

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

    echo "apptainer build oscar-count.sif recipe_oscar_count.sif"

    apptainer build oscar-count.sif recipe_oscar_count.sif

    echo ""

elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
fi

echo "Do you want to build the QC file? (Y/N)"
read -r choice

while [[ ! $choice =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter Y or N."
    read -r choice
done

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

    echo ""

    echo "apptainer build oscar-qc.sif recipe_oscar_qc.sif"

    apptainer build oscar-qc.sif recipe_oscar_qc.sif

elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
fi
