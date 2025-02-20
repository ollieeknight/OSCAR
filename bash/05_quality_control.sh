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

IFS=',' read -r -a project_ids <<< "${project_id}"

output_project_id="${project_ids[0]}"
output_project_dir="${dir_prefix}/${output_project_id}"
output_project_scripts="${output_project_dir}/${output_project_id}_scripts"
output_project_libraries="${output_project_scripts}/libraries"
output_project_outs="${output_project_dir}/${output_project_id}_outs"

# Check if metadata file exists for all project IDs
for project_id in "${project_ids[@]}"; do
    project_dir="${dir_prefix}/${project_id}"
    project_scripts="${project_dir}/${project_id}_scripts"
    metadata_file="${project_scripts}/metadata/${metadata_file_name}"
    
    if [ ! -f "${metadata_file}" ]; then
        echo -e "\033[0;31mERROR:\033[0m Metadata file for ${project_id} not found, please check path"
        exit 1
    fi
done

check_folder_exists "${output_project_outs}"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers

qc_container=${TMPDIR}/OSCAR/oscar-qc_latest.sif

metadata_file="${output_project_scripts}/metadata/metadata.csv"

libraries=($(find ${output+project_outs}/ -maxdepth 1 -mindepth 1 -type d -not -name 'logs' -exec basename {} \;))
mapfile -t libraries < <(printf '%s\n' "${libraries[@]}" | sort)

for library in "${libraries[@]}"; do

        read assay experiment_id historical_number replicate modality < <(extract_variables "$library")

        n_donors=$(extract_donor_number_from_all_metadata "$library")

        feature_matrix_path=$(find "${output_project_outs}/${library}/" -type f -name "raw_feature_bc_matrix.h5" -print -quit)
        peak_matrix_path=$(find "${output_project_outs}/${library}/" -type f -name "raw_peak_bc_matrix.h5" -print -quit)

        if [ -n "$feature_matrix_path" ]; then
                read -p "Submit ambient RNA removal with cellbender for ${library}? (Y/N): " choice
                while [[ ! $choice =~ ^[YyNn]$ ]]; do
                        echo "Invalid input. Please enter Y or N."
                        read -p "Submit ambient RNA removal with cellbender for ${library}? (Y/N): " choice
                done
                # Process choices
                if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
# Submit the job and capture job ID
job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name cellbender_${experiment_id}
#SBATCH --output ${output_project_outs}/logs/cellbender_${library}.out
#SBATCH --error ${output_project_outs}/logs/cellbender_${library}.out
#SBATCH --ntasks 1
#SBATCH --partition "gpu"
#SBATCH --gres gpu:1
#SBATCH --cpus-per-task 16
#SBATCH --mem 96GB
#SBATCH --time 18:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

log "OSCAR step 5: Quality control with cellbender"
log "See https://github.com/ollieeknight/OSCAR for more information"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Feature matrix          | ${feature_matrix_path}"
log "Output path             | ${output_project_outs}/${library}/cellbender"
log "Output file             | output.h5"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}

mkdir -p ${output_project_outs}/${library}/cellbender

echo ""

# Run cellbender
log "Starting CellBender remove-background..."
apptainer run --nv -B /data ${qc_container} cellbender remove-background \
        --cuda \
        --input ${feature_matrix_path} \
        --output ${output_project_outs}/${library}/cellbender/output.h5

echo ""

# Cleanup
rm ckpt.tar.gz

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
                        echo "Number of donors is $n_donors"
                        read -p "Would you like to genotype ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to genotype ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
# Submit the dependent job
sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${output_project_outs}/logs/geno_${library}.out
#SBATCH --error ${output_project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

log "OSCAR step 5: Quality control; genotyping with cellsnp-lite and vireo"
log "See https://github.com/ollieeknight/OSCAR for more information"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Input sample          | ${output_project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam"
log "Input cell barcodes   | ${output_project_outs}/${library}/cellbender/output_cell_barcodes.csv"
log "Output folder         | ${output_project_outs}/${library}/vireo"
log "VCF file              | /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
log "Number of donors      | $n_donors"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}
mkdir -p ${output_project_outs}/${library}/vireo

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
apptainer exec -B /data ${qc_container} cellsnp-lite \
        -s ${output_project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
        -b ${output_project_outs}/${library}/cellbender/output_cell_barcodes.csv \
        -O ${output_project_outs}/${library}/vireo \
        -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
        --minMAF 0.1 \
        --minCOUNT 20 \
        --gzip \
        -p \$(nproc)

echo ""

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${output_project_outs}/${library}/vireo \
        -o ${output_project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)

echo ""

log "All processing completed successfully!"
EOF
                                job_id=""
                        else
                                echo "Skipping genotyping for ${library}"
                        fi
                elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" == "" ]]; then
                        echo "Number of donors is $n_donors"
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
#SBATCH --output ${output_project_outs}/logs/geno_${library}.out
#SBATCH --error ${output_project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=32GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Input sample          | ${output_project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam"
log "Input cell barcodes   | ${output_project_outs}/${library}/cellbender/output_cell_barcodes.csv"
log "Output folder         | ${output_project_outs}/${library}/vireo"
log "VCF file              | /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
log "Number of donors      | $n_donors"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}
mkdir -p ${output_project_outs}/${library}/vireo

log ""

# Run cellsnp-lite
log "Starting cellsnp-lite processing..."
# apptainer exec -B /data ${qc_container} cellsnp-lite \
#         -s ${output_project_outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam \
#         -b ${output_project_outs}/${library}/cellbender/output_cell_barcodes.csv \
#         -O ${output_project_outs}/${library}/vireo \
#         -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
#         --minMAF 0.1 \
#         --minCOUNT 20 \
#         --gzip \
#         -p \$(nproc)

log ""

# Run vireo
log "Starting vireo processing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${output_project_outs}/${library}/vireo \
        -o ${output_project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)

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
                        echo "Number of donors is $n_donors"

                        read -p "Would you like to genotype ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to genotype ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${output_project_outs}/logs/geno_${library}.out
#SBATCH --error ${output_project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=128GB
#SBATCH --time=96:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Input sample          | ${output_project_outs}/${library}/outs/possorted_bam.bam"
log "Output name           | output"
log "Output folder         | ${output_project_outs}/${library}/mgatk"
log "Input barcodes             | ${output_project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv"
log "Number of donors      | $n_donors"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}

echo ""

# Run mgatk mtDNA genotyping
log "Starting mgatk mtDNA genotyping..."
# apptainer exec -B /data,/usr ${qc_container} mgatk tenx \
#         -i ${output_project_outs}/${library}/outs/possorted_bam.bam \
#         -n output \
#         -o ${output_project_outs}/${library}/mgatk \
#         -c 1 \
#         -bt CB \
#         -b ${output_project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
#         --skip-R

echo ""

rm -r ${output_project_outs}/${library}/.snakemake

mkdir -p ${output_project_outs}/${library}/AMULET

# Run AMULET doublet detection
log "Starting AMULET doublet detection..."
apptainer run -B /data ${qc_container} AMULET \
        ${output_project_outs}/${library}/outs/fragments.tsv.gz \
        ${output_project_outs}/${library}/outs/singlecell.csv \
        /opt/AMULET/human_autosomes.txt \
        /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed \
        ${output_project_outs}/${library}/AMULET \
        /opt/AMULET/

echo ""

mkdir -p ${output_project_outs}/${library}/vireo

# Run cellsnp-lite
log "Starting donor SNP genotyping with cellsnp-lite..."
# apptainer exec -B /data ${qc_container} cellsnp-lite \
#         -s ${output_project_outs}/${library}/outs/possorted_bam.bam \
#         -b ${output_project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
#         -O ${output_project_outs}/${library}/vireo \
#         -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
#         --minMAF 0.1 \
#         --minCOUNT 20 \
#         --gzip \
#         -p \$(nproc) \
#         --UMItag None

echo ""

# Run vireo
log "Starting vireo donor demultiplexing..."
apptainer run -B /data ${qc_container} vireo \
        -c ${output_project_outs}/${library}/vireo \
        -o ${output_project_outs}/${library}/vireo \
        -N $n_donors \
        -p \$(nproc)

EOF
                        else
                                echo "Skipping genotyping"
                        fi
                elif [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
                        echo "Number of donors is $n_donors"
                        read -p "Would you like to perform mitochondrial genotyping for ${library}? (Y/N): " choice
                        while [[ ! $choice =~ ^[YyNn]$ ]]; do
                                echo "Invalid input. Please enter Y or N."
                                read -p "Would you like to perform mitochondrial genotyping for ${library}? (Y/N): " choice
                        done
                        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name geno_${experiment_id}
#SBATCH --output ${output_project_outs}/logs/geno_${library}.out
#SBATCH --error ${output_project_outs}/logs/geno_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=128GB
#SBATCH --time=48:00:00

# Source the functions
source "${oscar_dir}/functions.sh"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Input sample          | ${output_project_outs}/${library}/outs/possorted_bam.bam"
log "Output name           | output"
log "Output folder         | ${output_project_outs}/${library}/mgatk"
log "Input barcodes             | ${output_project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv"
log "Number of donors      | $n_donors"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}

# Run mgatk mtDNA genotyping
log "Starting mgatk mtDNA genotyping..."
# apptainer exec -B /data,/usr ${qc_container} mgatk tenx \
#         -i ${output_project_outs}/${library}/outs/possorted_bam.bam \
#         -n output \
#         -o ${output_project_outs}/${library}/mgatk \
#         -c 8 \
#         -bt CB \
#         -b ${output_project_outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv \
#         --skip-R

echo ""

rm -r ${output_project_outs}/${library}/.snakemake

mkdir -p ${output_project_outs}/${library}/AMULET

# Run AMULET doublet detection
log "Starting AMULET doublet detection..."
apptainer run -B /data ${qc_container} AMULET \
        ${output_project_outs}/${library}/outs/fragments.tsv.gz \
        ${output_project_outs}/${library}/outs/singlecell.csv \
        /opt/AMULET/human_autosomes.txt \
        /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed \
        ${output_project_outs}/${library}/AMULET \
        /opt/AMULET/

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