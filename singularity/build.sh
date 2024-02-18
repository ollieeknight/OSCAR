#!/bin/bash/

if [[ ! -f ../cellranger-7.2.0.tar.gz && ! -f ../cellranger-atac-2.1.0.tar.gz ]]; then
    echo "First, you need to download cellranger-7.2.0.tar.gz AND cellranger-atac-2.1.0.tar.gz ABOVE this folder, then come back."
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

    cd ..

    echo "apptainer build oscar-count.sif singularity/recipe_oscar-count.sif"

    apptainer build singularity/oscar-count.sif singularity/recipe_oscar-count.sif

    echo ""

    cd singularity

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

    cd ..

    echo ""

    echo "apptainer build oscar-qc.sif singularity/recipe_oscar-qc.sif"

    apptainer build singularity/oscar-qc.sif singularity/recipe_oscar-qc.sif

    cd singularity

elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
fi
