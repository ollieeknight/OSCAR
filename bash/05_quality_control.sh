#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"

# Function to display help message
display_help() {
    echo "Usage: ${0} [options]"
    echo ""
    echo "Options:"
    echo "  --project-id <id>                             Set the project ID. Can be multiple (comma-separated)"
    echo "  --dir-prefix <path>                           Set the directory prefix (default: ${HOME}/scratch/ngs)"
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
project_scripts="${project_dir}/${project_id}_scripts"

check_folder_exists "${project_outs}"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers
qc_container=${TMPDIR}/OSCAR/oscar-qc_latest.sif

metadata_file="${project_dir}/${project_id}_scripts/metadata/metadata.csv"

libraries=($(find ${project_outs}/ -maxdepth 1 -mindepth 1 -type d -not -name 'logs' -exec basename {} \;))
mapfile -t libraries < <(printf '%s\n' "${libraries[@]}" | sort)

for library in "${libraries[@]}"; do

        read assay experiment_id historical_number replicate modality < <(extract_variables "$library")
        log "library: ${library}"
        log "assay: ${assay}"
        log "experiment_id: ${experiment_id}"
        log "historical_number: ${historical_number}"
        log "replicate: ${replicate}"

        read n_donors < <(search_metadata "$library" "$assay" "$experiment_id" "$historical_number" "$replicate")
        log "n_donors: ${n_donors}"
        feature_matrix_path=$(find "${project_outs}/${library}/" -type f -name "raw_feature_bc_matrix.h5" -print -quit)
        peak_matrix_path=$(find "${project_outs}/${library}/" -type f -name "raw_peak_bc_matrix.h5" -print -quit)

        if [ -n "$feature_matrix_path" ]; then
                read -p "Submit ambient RNA removal with cellbender for ${library}? (Y/N): " choice
                while [[ ! $choice =~ ^[YyNn]$ ]]; do
                        echo "Invalid input. Please enter Y or N."
                        read -p "Submit ambient RNA removal with cellbender for ${library}? (Y/N): " choice
                done
                # Process choices
                if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
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

# Log input variables
log "Input variables:"
log "experiment_id: ${experiment_id}"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "qc_container: ${qc_container}"
log "feature_matrix_path: ${feature_matrix_path}"

log ""

# Change to library directory
log "Changing to library directory..."
cd ${project_outs}/${library}
check_status "Directory change"

log ""

# Create cellbender directory
log "Creating cellbender directory..."
mkdir -p ${project_outs}/${library}/cellbender
check_status "Directory creation"

log ""

# Run cellbender
log "Starting CellBender remove-background..."
apptainer run --nv -B /data ${qc_container} cellbender remove-background \
        --cuda \
        --input ${feature_matrix_path} \
        --output ${project_outs}/${library}/cellbender/output.h5
check_status "CellBender processing"

log ""

# Cleanup
log "Cleaning up temporary files..."
rm ckpt.tar.gz
check_status "Cleanup"

log ""

log "CellBender processing completed successfully!"
EOF
)
                job_id=$(echo "$job_id" | awk '{print $4}')
                elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                        continue
                fi

                if [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
                        echo "Skipping genotyping for ${library}, as this is either a mouse run, or only contains 1 donor"
                        job_id=""
                elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" != "" ]]; then
                        read -p "Would you like to genotype ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to genotype ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
                                echo "Submitting vireo genotyping for ${library}"
# Submit the dependent job
sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${project_outs}/logs/geno_${library}.out
#SBATCH --error ${project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

# Log input variables
log "Input variables:"
log "experiment_id: ${experiment_id}"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "qc_container: ${qc_container}"
log "n_donors: ${n_donors}"

log ""

cd ${project_outs}/${library}
mkdir -p ${project_outs}/${library}/vireo

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
apptainer exec -B /data ${qc_container} cellsnp-lite \
        -s ${project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
        -b ${project_outs}/${library}/cellbender/output_cell_barcodes.csv \
        -O ${project_outs}/${library}/vireo \
        -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip \
        -p \$(nproc)
check_status "Cellsnp-lite processing"

log ""

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${project_outs}/${library}/vireo \
        -o ${project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)
check_status "Vireo processing"

log ""

log "All processing completed successfully!"
EOF
                                job_id=""
                        else
                                echo "Skipping genotyping for ${library}"
                        fi
                elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" == "" ]]; then
                        read -p "Would you like to genotype ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to genotype ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
# Submit the job
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${project_outs}/logs/geno_${library}.out
#SBATCH --error ${project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

# Log input variables
log "Input variables:"
log "experiment_id: ${experiment_id}"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "qc_container: ${qc_container}"
log "n_donors: ${n_donors}"

log ""

cd ${project_outs}/${library}
mkdir -p ${project_outs}/${library}/vireo

log ""

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
apptainer exec -B /data ${qc_container} cellsnp-lite \
        -s ${project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
        -b ${project_outs}/${library}/cellbender/output_cell_barcodes.csv \
        -O ${project_outs}/${library}/vireo \
        -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip \
        -p \$(nproc)
check_status "Cellsnp-lite processing"

log ""

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${project_outs}/${library}/vireo \
        -o ${project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)
check_status "Vireo processing"

log ""

log "All processing completed successfully!"
EOF
                                job_id=""
                        else
                                echo "Skipping genotyping for ${library}"
                        fi
                else
                        echo -e "\033[0;31mERROR:\033[0m Cannot determine the number of donors for ${library}"
                        exit 1
                fi
        elif [ -n "$peak_matrix_path" ]; then
                if [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' ]]; then
                        read -p "Would you like to genotype ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to genotype ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
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

cd ${project_outs}/${library}

log ""

# Run mgatk mtDNA genotyping
log "Starting mgatk mtDNA genotyping..."
apptainer exec -B /data,/usr ${qc_container} mgatk tenx \
        -i ${project_outs}/${library}/outs/possorted_bam.bam \
        -n output \
        -o ${project_outs}/${library}/mgatk \
        -c 1 \
        -bt CB \
        -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
        --skip-R
check_status "mgatk processing"

log ""

rm -r ${project_outs}/${library}/.snakemake

mkdir -p ${project_outs}/${library}/AMULET

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

log ""

mkdir -p ${project_outs}/${library}/vireo

# Run cellsnp-lite
log "Starting donor SNP genotyping with cellsnp-lite..."
apptainer exec -B /data ${qc_container} cellsnp-lite \
        -s ${project_outs}/${library}/outs/possorted_bam.bam \
        -b ${project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
        -O ${project_outs}/${library}/vireo \
        -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip \
        -p \$(nproc) \
        --UMItag None
check_status "Cellsnp-lite processing"

log ""

# Run vireo
log "Starting vireo donor demultiplexing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${project_outs}/${library}/vireo \
        -o ${project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)
check_status "Vireo processing"

log ""

log "All processing completed successfully!"
EOF
                        else
                                echo "Skipping genotyping"
                        fi
                elif [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
                        read -p "Would you like to perform mitochondrial genotyping for ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to perform mitochondrial genotyping for ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
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

# Log input variables
log "Input variables:"
log "experiment_id: ${experiment_id}"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "qc_container: ${qc_container}"
log "n_donors: ${n_donors}"

log ""

cd ${project_outs}/${library}

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

log ""

rm -r ${project_outs}/${library}/.snakemake

mkdir -p ${project_outs}/${library}/AMULET

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

log ""

log "Processing completed successfully!"
EOF
                        else
                                echo "Skipping genotyping"
                        fi
                fi
        else
                # Action when neither file is found
                echo -e "\033[0;31mERROR:\033[0m Neither feature matrix nor peak matrix was found for ${library}"
                exit 1
        fi
done