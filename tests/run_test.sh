#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --reference <reference_path>"
    echo "  --reference    Path to the human reference transcriptome"
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --reference)
            reference="$2"
            shift 2
            ;;
        *)
            echo "Invalid option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check if reference is provided
if [ -z "$reference" ]; then
    echo "Error: Please provide a human reference transcriptome with the --reference option."
    usage
    exit 1
fi

num_cores=$(nproc)

echo ""
echo "Downloading required run files..."
echo ""

# Create directories and download required files
mkdir -p test_run/
wget -q -O test_run/cellranger-tiny-bcl-1.2.0.tar.gz \
    https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz

tar -xf test_run/cellranger-tiny-bcl-1.2.0.tar.gz -C test_run/
mv test_run/cellranger-tiny-bcl-1.2.0 test_run/test_run_bcl
rm test_run/cellranger-tiny-bcl-1.2.0.tar.gz

mkdir -p test_run/test_run_scripts/indices/
wget -q -O test_run/test_run_scripts/indices/test_run.csv \
    https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv

# Pull OSCAR count image
oscar_image_path="${TMPDIR}/OSCAR/oscar-count_latest.sif"
if [ ! -f "$oscar_image_path" ]; then
    echo "Pulling OSCAR count image to ${TMPDIR}/OSCAR, this might take some time..."
    mkdir -p "${TMPDIR}/OSCAR"
    apptainer pull "$oscar_image_path" library://romagnanilab/oscar/oscar-count:latest
fi

cd test_run || exit
echo "Running cellranger mkfastq..."
echo ""
apptainer run -B /data "$oscar_image_path" \
    cellranger mkfastq --id test_run_fastq --run test_run_bcl/ --csv test_run_scripts/indices/test_run.csv &> mkfastq.log

echo "Running cellranger count..."
echo ""
apptainer run -B /data "$oscar_image_path" \
    cellranger count --id test_run_sample --fastqs test_run_fastq/outs/fastq_path/H35KCBCXY/ \
    --sample test_sample --localcores "$num_cores" --transcriptome "$reference" \
    --chemistry SC3Pv3 --no-bam &> count.log

# Check if output directory exists
if [ -d "test_run_sample/outs/" ]; then
    echo "Success! Shutting down."
else
    echo "Uh oh! Output directory not found."
fi
