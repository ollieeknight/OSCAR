#! /bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --reference <reference_path>"
    echo "  --reference    Path to the human reference transcriptome"
}

# Parse command line arguments
while [[ "$@" ]]; do
    case "$1" in
        --reference)
            if [ -n "$2" ]; then
                reference="$2"
            else
                echo "Error: --reference option requires a non-empty argument."
            echo "Invalid option: $1" >&2
            exit 1
            fi
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

echo -e "\nDownloading required run files...\n"
echo ""

# Create directories and download required files
mkdir -p test_run/
wget -q -O test_run/cellranger-tiny-bcl-1.2.0.tar.gz \
    https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
if [ ! -f "test_run/cellranger-tiny-bcl-1.2.0.tar.gz" ]; then
    wget -q -O test_run/cellranger-tiny-bcl-1.2.0.tar.gz \
        https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
fi
mv test_run/cellranger-tiny-bcl-1.2.0 test_run/test_run_bcl
rm test_run/cellranger-tiny-bcl-1.2.0.tar.gz

mkdir -p test_run/test_run_scripts/indices/
wget -q -O test_run/test_run_scripts/indices/test_run.csv \
    https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
oscar_count_image_path="${TMPDIR}/OSCAR/oscar-count_latest.sif"
if [ ! -f "$oscar_count_image_path" ]; then
    echo "Pulling OSCAR count image to ${TMPDIR}/OSCAR, this might take some time..."
    mkdir -p "${TMPDIR}/OSCAR"
    apptainer pull "$oscar_count_image_path" library://romagnanilab/oscar/oscar-count:latest
apptainer pull "$oscar_image_path" library://romagnanilab/oscar/oscar-count:latest
fi

# Run the OSCAR count image
apptainer run -B /data "$oscar_count_image_path"

echo ""

# Run cellranger mkfastq to generate FASTQ files from BCL files
echo "Running cellranger count..."
if [ ! -f "$oscar_image_path" ]; then
    echo "Error: OSCAR image not found at $oscar_image_path"
    exit 1
fi
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
