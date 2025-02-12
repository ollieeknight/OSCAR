#!/bin/bash

# Set the directory of the current script
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
# Source the functions script
source "${oscar_dir}/functions.sh"
# Set default directory prefix and metadata file name
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

# Function to display help message
display_help() {
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  --project-id <id>          Set the project ID"
  echo "  --dir-prefix <path>        Set the directory prefix (default: ${HOME}/scratch/ngs)"
  echo "  --metadata-file-name <name> Set the metadata file name (default: metadata.csv)"
  echo "  --help                     Display this help message"
  exit 0
}

# Parse command line arguments
while [[ "${#}" -gt 0 ]]; do
  if [[ "${1}" == --* ]]; then
    if [[ "${1}" == "--help" ]]; then
      display_help  # Display help message and exit
    fi
    if [[ -z "${2}" ]]; then
      echo "Error: Missing value for parameter ${1}"
      exit 1
    fi
    var_name=$(echo "${1}" | sed 's/--//; s/-/_/')
    declare "${var_name}"="${2}"
    shift 2
  else
    echo "Invalid option: ${1}"
    exit 1
  fi
done

check_project_id

# Set project directories
project_dir="${dir_prefix}/${project_id}"
project_scripts="${project_dir}/${project_id}_scripts"
project_indices="${project_scripts}/indices"

# Check if metadata file exists
metadata_file="${project_scripts}/metadata/${metadata_file_name}"
check_metadata_file "${metadata_file}"

# Determine the run type
run_type=$(check_run_type "${project_id}" "${dir_prefix}")

# Check if indices folder exists
check_folder_exists "${project_scripts}/indices"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers

count_container=$TMPDIR/OSCAR/oscar-count_latest.sif

# Check base masks for step 1 and step 2
check_base_masks_step1
check_base_masks_step2

# List index files and extract flowcell ID from RunInfo.xml
index_files=($(ls "${project_dir}/${project_id}_scripts/indices"))
flowcell_id=$(grep "<Flowcell>" "${project_dir}/${project_id}_bcl/RunInfo.xml" | sed -e 's|.*<Flowcell>\(.*\)</Flowcell>.*|\1|')

# Loop through each index file
for file in "${index_files[@]}"; do
    index_file="${file%.*}"

    # Check base masks for step 3 and parse the command, index type, filter option, and base mask
    read -r cellranger_command index_type filter_option base_mask < <(check_base_masks_step3 "$file" "$run_type")
    cellranger_command="${cellranger_command//./ }"
    index_type="${index_type//./ }"
    filter_option="${filter_option//./ }"
    base_mask="${base_mask//./ }"

    # Prompt the user for confirmation
    read -p $'\033[0;33mINPUT REQUIRED:\033[0m Submit '"${index_file}"' for fastq demultiplexing? (Y/N) ' choice
    while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
        echo "Invalid input. Please enter y or n"
        read -p $'\033[0;33mINPUT REQUIRED:\033[0m Submit '"${index_file}"' for fastq demultiplexing? (Y/N) ' choice
    done

    # If user confirms, submit the job to SLURM
    if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

        # Create logs directory
        mkdir -p "${project_dir}/${project_id}_fastq/logs/"

        # Submit the job to SLURM
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${project_id}
#SBATCH --output ${project_dir}/${project_id}_fastq/logs/${index_file}.out
#SBATCH --error ${project_dir}/${project_id}_fastq/logs/${index_file}.out
#SBATCH --ntasks=16
#SBATCH --mem=32GB
#SBATCH --time=4:00:00

source "${oscar_dir}/functions.sh"

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Cellranger flavour      | ${cellranger_command}"
log "Sequencing run          | ${project_dir}/${project_id}_bcl"
log "Flowcell ID             | ${flowcell_id}"
log "Index .csv file  | ${project_scripts}/indices/${file}"
log "Output name             | ${index_file}"
log "Base mask               | ${base_mask}"
log "Filter option           | ${filter_option}"
log "----------------------------------------"

mkdir -p ${project_dir}/${project_id}_fastq/logs
cd ${project_dir}/${project_id}_fastq/

echo ""

# Run cellranger
log "Running ${cellranger_command}"

echo ""

apptainer run -B /data ${count_container} \
    ${cellranger_command} \
    --run ${project_dir}/${project_id}_bcl \
    --id ${index_file} \
    --csv ${project_scripts}/indices/${file} \
    --use-bases-mask ${base_mask} \
    --delete-undetermined \
    --barcode-mismatches 1 \
    ${filter_option}

log ""

mkdir -p ${project_dir}/${project_id}_fastq/${index_file}/fastqc

# Run fastqc
log "Running FastQC"
find "${project_dir}/${project_id}_fastq/${index_file}/outs/fastq_path/${flowcell_id}"* -name "*.fastq.gz" | \
    parallel -j \$(nproc) "apptainer run -B /data ${count_container} fastqc {} \
    --outdir ${project_dir}/${project_id}_fastq/${index_file}/fastqc"

echo ""

check_status "fastqc"

echo ""

# Run multiqc to aggregate fastqc reports
log "Running MultiQC"
apptainer run -B /data ${count_container} multiqc \
    "${project_dir}/${project_id}_fastq/${index_file}" \
    -o "${project_dir}/${project_id}_fastq/${index_file}/multiqc"

echo ""

# Clean up temporary files
log "Cleaning up temporary files"
rm -r ${project_dir}/${project_id}_fastq/${index_file}/_* ${project_dir}/${project_id}_fastq/${index_file}/MAKE*
echo ""

log "All processing completed successfully"
EOF
    elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
        :
    else
        echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
    fi
done
