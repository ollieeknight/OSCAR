#!/bin/bash/

echo "First, you need to download cellranger-7.2.0.tar.gz AND cellranger-atac-2.1.0.tar.gz ABOVE this folder, then come back."
echo "You can download them from here:"
echo "https://www.10xgenomics.com/support/software/cell-ranger/downloads"
echo "https://support.10xgenomics.com/single-cell-atac/software/downloads/latest"

echo "Are you sure you want to build these singularity files? (Y/N)"
read -r choice

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

	cd ..

	echo "apptainer build oscar-counting.sif singularity/recipe_oscar-counting.sif"
	apptainer build oscar-counting.sif singularity/recipe_oscar-counting.sif
	echo ""

        echo "apptainer build oscar-qc.sif singularity/recipe_oscar-qc.sif"
        apptainer build oscar-qc.sif singularity/recipe_oscar-qc.sif
elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
else
    echo "Invalid choice. No worries"
fi
