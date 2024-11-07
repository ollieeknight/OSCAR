# Set the directory of the current script
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
# Source the functions script
source "${oscar_dir}/functions.sh"
# Set default directory prefix and metadata file name
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

# Check if project_id is set
check_project_id

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
  if [[ "$1" == --* ]]; then
    var_name=$(echo "$1" | sed 's/--//; s/-/_/')
    declare "$var_name"="$2"
    shift 2
  else
    echo "Invalid option: $1"
    exit 1
  fi
done

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

    # Print the command to be executed
    echo ""
    echo "apptainer run -B /data ${container} ${cellranger_command} --id ${index_file} --run ${project_dir}/${project_id}_bcl --csv ${project_scripts}/indices/${file} --use-bases-mask ${base_mask} --delete-undetermined --barcode-mismatches 1 ${filter_option}"
    echo ""

    # Prompt the user for confirmation
    echo -e "\033[0;33mINPUT REQUIRED:\033[0m Is this alright? (Y/N)"
    read -r choice
    while [[ ! $choice =~ ^[YyNn]$ ]]; do
        echo "Invalid input. Please enter y or n"
        read -r choice
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
#SBATCH --mem=64000
#SBATCH --time=12:00:00

# Change to the fastq directory
cd ${project_dir}/${project_id}_fastq/

# Run the cellranger command
apptainer run -B /data ${container} ${cellranger_command} --id ${index_file} --run ${project_dir}/${project_id}_bcl --csv ${project_scripts}/indices/${file} --use-bases-mask ${base_mask} --delete-undetermined --barcode-mismatches 1 ${filter_option}

# Create fastqc directory
mkdir -p ${project_dir}/${project_id}_fastq/${index_file}/fastqc

# Run fastqc on all fastq.gz files in parallel
find "${project_dir}/${project_id}_fastq/${index_file}/outs/fastq_path/${flowcell_id}"* -name "*.fastq.gz" | parallel -j $(nproc) "apptainer run -B /data ${container} fastqc {} --outdir ${project_dir}/${project_id}_fastq/${index_file}/fastqc"

# Run multiqc to aggregate fastqc reports
apptainer run -B /data ${container} multiqc "${project_dir}/${project_id}_fastq/${index_file}" -o "${project_dir}/${project_id}_fastq/${index_file}/multiqc"

# Clean up temporary files
rm -r ${project_dir}/${project_id}_fastq/${index_file}/_* ${project_dir}/${project_id}_fastq/${index_file}/MAKE*
EOF
    elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
        :
    else
        echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
    fi
done
