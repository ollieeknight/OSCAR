#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

# Function to display help message
display_help() {
  echo "Usage: ${0} [options]"
  echo ""
  echo "Options:"
  echo "  --project-id <id>                             Set the project ID. Can be multiple (comma-separated)"
  echo "  --dir-prefix <path>                           Set the directory prefix (default: ${HOME}/scratch/ngs)"
  echo "  --metadata-file-name <name>                   Set the metadata file name (default: metadata.csv)"
  echo "  --help                                        Display this help message"
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

# Define project directories
project_dir="${dir_prefix}/${project_id}"
project_outs="${project_dir}/${project_id}_outs"

check_folder_exists "${project_outs}"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers
qc_container=${TMPDIR}/OSCAR/oscar-qc_latest.sif

metadata_file="${project_dir}/${project_id}_scripts/metadata/metadata.csv"

libraries=($(find ${project_outs}/ -maxdepth 1 -mindepth 1 -type d -not -name 'logs' -exec basename {} \;))
mapfile -t libraries < <(printf '%s\n' "${libraries[@]}" | sort)

for library in "${libraries[@]}"; do

    read assay experiment_id historical_number replicate < <(extract_variables "$library")

    read n_donors ADT_file < <(search_metadata "$library" "$assay" "$experiment_id" "$historical_number" "$replicate" project_ids[@] "$dir_prefix")

    feature_matrix_path=$(find "${project_outs}/${library}/" -type f -name "raw_feature_bc_matrix.h5" -print -quit)
    peak_matrix_path=$(find "${project_outs}/${library}/" -type f -name "raw_peak_bc_matrix.h5" -print -quit)

    if [ -n "$feature_matrix_path" ]; then
        read -p "Would you like to submit ambient RNA removal with cellbender for ${library}? (Y/N)" perform_function
        # Convert input to uppercase for case-insensitive comparison
        perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
        if [ "$perform_function" != "Y" ]; then
            echo "Skipping cellbender for ${library}"
        else
                echo "Submitting cellbender for ${library}"
# Submit the job and capture job ID
job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name cellbender_${experiment_id}
#SBATCH --output ${project_outs}/logs/cellbender_${library}.out
#SBATCH --error ${project_outs}/logs/cellbender_${library}.out
#SBATCH --ntasks 1
#SBATCH --partition "gpu"
#SBATCH --gres gpu:1
#SBATCH --cpus-per-task 16
#SBATCH --mem 96GB
#SBATCH --time 18:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

# Create cellbender directory
log "Creating cellbender directory..."
mkdir -p ${project_outs}/${library}/cellbender
check_status "Directory creation"

# Run cellbender
log "Starting CellBender remove-background..."
apptainer run --nv -B /data ${qc_container} cellbender remove-background \
    --cuda \
    --input ${feature_matrix_path} \
    --output ${project_outs}/${library}/cellbender/output.h5
check_status "CellBender processing"

# Cleanup
log "Cleaning up temporary files..."
rm ckpt.tar.gz
check_status "Cleanup"

log "CellBender processing completed successfully!"
EOF
)
        )
        job_id=$(echo "$job_id" | awk '{print $4}')
        echo ""
        fi
        if [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            echo "Skipping genotyping for ${library}, as this is either a mouse run, or only contains 1 donor"
            job_id=""
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" != "" ]]; then
             read -p "Would you like to genotype ${library}? (Y/N)" perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping for ${library}"
            else
                echo "Submitting vireo genotyping for ${library}"
# Submit the dependent job
sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name vireo_${experiment_id}
#SBATCH --output ${project_outs}/logs/vireo_${library}.out
#SBATCH --error ${project_outs}/logs/vireo_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00
# The following line ensures that this job runs after the previous job with ID $job_id
#SBATCH --dependency=afterok:$job_id

# Source the functions
source "${oscar_dir}/functions.sh"

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

# Create vireo directory
log "Creating vireo directory..."
mkdir -p ${project_outs}/${library}/vireo
check_status "Directory creation"

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
apptainer run -B /data ${qc_container} cellsnp-lite \
    -s ${project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
    -b ${project_outs}/${library}/cellbender/output_cell_barcodes.csv \
    -O ${project_outs}/${library}/vireo \
    -R /data/cephfs-2/unmirrored/groups/romagnani/work/ref/vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    --minMAF 0.1 \
    --minCOUNT 20 \
    --gzip \
    -p \$(nproc)
check_status "Cellsnp-lite processing"

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
    -c ${project_outs}/${library}/vireo \
    -o ${project_outs}/${library}/vireo \
    -N $n_donors \
    -p \$(nproc)
check_status "Vireo processing"

log "All processing completed successfully!"
EOF
                job_id=""
            fi
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" == "" ]]; then
             read -p "Would you like to genotype ${library}? (Y/N)" perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping for ${library}"
            else
                echo "Submitting vireo genotyping for ${library}"
# Submit the job
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name vireo_${experiment_id}
#SBATCH --output ${project_outs}/logs/vireo_${library}.out
#SBATCH --error ${project_outs}/logs/vireo_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

# Create vireo directory
log "Creating vireo directory..."
mkdir -p ${project_outs}/${library}/vireo
check_status "Directory creation"

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
apptainer run -B /data ${qc_container} cellsnp-lite \
    -s ${project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
    -b ${project_outs}/${library}/cellbender/output_cell_barcodes.csv \
    -O ${project_outs}/${library}/vireo \
    -R /data/cephfs-2/unmirrored/groups/romagnani/work/ref/vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    --minMAF 0.1 \
    --minCOUNT 20 \
    --gzip \
    -p \$(nproc)
check_status "Cellsnp-lite processing"

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
    -c ${project_outs}/${library}/vireo \
    -o ${project_outs}/${library}/vireo \
    -N $n_donors \
    -p \$(nproc)
check_status "Vireo processing"

log "All processing completed successfully!"
EOF
                job_id=""
            fi
        else
            echo -e "\033[0;31mERROR:\033[0m Cannot determine the number of donors for ${library}"
            exit 1
        fi
    elif [ -n "$peak_matrix_path" ]; then
        if [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' ]]; then
            read -p "Would you like to genotype ${library}? (Y/N): " perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')

            # Check if the input is 'N' or 'n'
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping"
            else
                echo "Submitting vireo genotyping for ${library}"
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${project_outs}/logs/geno_${library}.out
#SBATCH --error ${project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=128GB
#SBATCH --time=96:00:00
# Source the functions
source "${oscar_dir}/functions.sh"

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

# Run mgatk mtDNA genotyping
log "Starting mgatk mtDNA genotyping..."
apptainer exec -B /data,/usr ${qc_container} mgatk tenx \
    -i ${project_outs}/${library}/outs/possorted_bam.bam \
    -n output \
    -o ${project_outs}/${library}/mgatk \
    -c 8 \
    -bt CB \
    -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
    --skip-R
check_status "mgatk processing"

# Cleanup snakemake
log "Cleaning up snakemake files..."
rm -r ${project_outs}/${library}/.snakemake
check_status "Snakemake cleanup"

# Create AMULET directory
log "Creating AMULET directory..."
mkdir -p ${project_outs}/${library}/AMULET
check_status "AMULET directory creation"

# Run AMULET doublet detection
log "Starting AMULET doublet detection..."
apptainer run -B /data ${qc_container} AMULET \
    ${project_outs}/${library}/outs/fragments.tsv.gz \
    ${project_outs}/${library}/outs/singlecell.csv \
    /opt/AMULET/human_autosomes.txt \
    /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed \
    ${project_outs}/${library}/AMULET \
    /opt/AMULET/
check_status "AMULET processing"

# Create vireo directory
log "Creating vireo directory..."
mkdir -p ${project_outs}/${library}/vireo
check_status "Vireo directory creation"

# Run cellsnp-lite
log "Starting donor SNP genotyping with cellsnp-lite..."
apptainer run -B /data ${qc_container} cellsnp-lite \
    -s ${project_outs}/${library}/outs/possorted_bam.bam \
    -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
    -O ${project_outs}/${library}/vireo \
    -R /data/cephfs-2/unmirrored/groups/romagnani/work/ref/vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    --minMAF 0.1 \
    --minCOUNT 20 \
    --gzip \
    -p \$(nproc) \
    --UMItag None
check_status "Cellsnp-lite processing"

# Run vireo
log "Starting vireo donor demultiplexing..."
apptainer run -B /data ${qc_container} vireo \
    -c ${project_outs}/${library}/vireo \
    -o ${project_outs}/${library}/vireo \
    -N $n_donors \
    -p \$(nproc)
check_status "Vireo processing"

log "All processing completed successfully!"
EOF
            fi
        elif [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            read -p "Would you like to perform mitochondrial genotyping for ${library}? (Y/N): " perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            # Check if the input is 'N' or 'n'
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping"
            else
                echo "Submitting genotyping for ${library}"
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${project_outs}/logs/geno_${library}.out
#SBATCH --error ${project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=32GB
#SBATCH --time=48:00:00
# Source the functions
source "${oscar_dir}/functions.sh"

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

# Run mgatk mtDNA genotyping
log "Starting mgatk mtDNA genotyping..."
apptainer exec -B /data,/usr ${qc_container} mgatk tenx \
    -i ${project_outs}/${library}/outs/possorted_bam.bam \
    -n output \
    -o ${project_outs}/${library}/mgatk \
    -c 8 \
    -bt CB \
    -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
    --skip-R
check_status "mgatk processing"

# Cleanup snakemake files
log "Cleaning up snakemake files..."
rm -r ${project_outs}/${library}/.snakemake
check_status "Snakemake cleanup"

# Create AMULET directory
log "Creating AMULET directory..."
mkdir -p ${project_outs}/${library}/AMULET
check_status "Directory creation"

# Run AMULET doublet detection
log "Starting AMULET doublet detection..."
apptainer run -B /data ${qc_container} AMULET \
    ${project_outs}/${library}/outs/fragments.tsv.gz \
    ${project_outs}/${library}/outs/singlecell.csv \
    /opt/AMULET/human_autosomes.txt \
    /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed \
    ${project_outs}/${library}/AMULET \
    /opt/AMULET/
check_status "AMULET processing"

log "Processing completed successfully!"
EOF
            fi
        fi
	else
        # Action when neither file is found
        echo -e "\033[0;31mERROR:\033[0m Neither feature matrix nor peak matrix was found for ${library}"
        exit 1
    fi
done