#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

# Function to display help message
display_help() {
    echo "Usage: $0 [options]"
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
project_scripts="${project_dir}/${project_id}_scripts"
project_indices="${project_scripts}/indices"
project_libraries="${project_scripts}/libraries"
project_outs="${project_dir}/${project_id}_outs"

# Check if metadata file exists
metadata_file="${project_scripts}/metadata/${metadata_file_name}"
check_metadata_file "${metadata_file}"

# Check if indices folder exists
check_folder_exists "${project_scripts}/indices"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers
count_container=${TMPDIR}/OSCAR/oscar-count_latest.sif

# Take the csv files into a list and remove the .csv suffix
libraries=($(ls "${project_libraries}" | awk -F/ '{print $NF}' | awk -F. '{print $1}'))
mkdir -p ${project_outs}/

# Iterate over each library file to submit counting jobs
for library in "${libraries[@]}"; do
    # Skip processing lines with 'ADT' in the library name
    if grep -q '.*\(HTO\|ADT\).*' "${project_libraries}/${library}.csv"; then
        echo "Processing ${library} as an ADT/HTO library"
        continue
    elif grep -q '.*\(ATAC\).*' "${project_libraries}/${library}.csv"; then
        echo "Processing ${library} as an ATAC library"
        fastq_names=""
        fastq_dirs=""

        read fastq_names fastq_dirs < <(count_read_csv "${project_libraries}" "$library")

        extra_arguments=$(count_check_dogma "${project_libraries}" "$library")

        read -p "Submit ${library} for cellranger-atac count? (Y/N): " choice
        while [[ ! $choice =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -p "Submit ${library} for cellranger-atac count? (Y/N): " choice
        done
        # Process choices
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
            mkdir -p ${project_outs}/logs
            # Submit the job to slurm for counting
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output ${project_outs}/logs/cellranger_${library}.out
#SBATCH --error ${project_outs}/logs/cellranger_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=64GB
#SBATCH --time=72:00:00

source "${oscar_dir}/functions.sh"

# Log input variables
log "Input variables:"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "count_container: ${count_container}"
log "fastq_dirs: ${fastq_dirs}"
log "fastq_names: ${fastq_names}"
log "extra_arguments: ${extra_arguments}"

log ""

cd ${project_outs}

# Run cellranger-atac count
log "Running cellranger-atac count"
apptainer run -B /data ${count_container} cellranger-atac count \
    --id $library \
    --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ \
    --fastqs $fastq_dirs \
    --sample $fastq_names \
    --localcores \$(nproc) \
    $extra_arguments
check_status "cellranger-atac count count"

rm -r ${project_outs}/$library/_* ${project_outs}/$library/SC_ATAC_COUNTER_CS

log "All processing completed successfully"
EOF
            )
            count_submitted='YES'
            job_id=$(echo "$job_id" | awk '{print $4}')
        elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
            count_submitted='NO'
        else
            echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
        fi

        if grep -q '.*\(ASAP\).*' "${project_libraries}/${library}.csv"; then

            # Extract the ADT file name from the metadata
            ADT_file=$(extract_adt_file "$metadata_file" "$library")

            # Check if the ADT file exists
            if [[ ! -f "${project_libraries}/${ADT_file}" ]]; then
                echo "ERROR: ${project_libraries}/${ADT_file} not found."
                exit 1
            fi
            echo "DEBUG: ${project_libraries}/${ADT_file} found."

            # Determine the correct ADT CSV file name by replacing _ATAC with _ADT
            adt_library_csv="${library/_ATAC/_ADT}.csv"

            echo "DEBUG: adt_library_csv is: ${adt_library_csv}"

            # Check if the ADT CSV file exists
            if [[ ! -f "${project_libraries}/${adt_library_csv}" ]]; then
                echo "ERROR: ${project_libraries}/${adt_library_csv} not found."
                exit 1
            fi
            echo "DEBUG: ${project_libraries}/${adt_library_csv} found."

            read -p "Perform ADT counting? (Y/N): " choice
            while [[ ! $choice =~ ^[YyNn]$ ]]; do
                echo "Invalid input. Please enter Y or N."
                read -p "Perform ADT counting? (Y/N): " choice
            done

            if [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                continue
            elif [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

                read fastq_dirs fastq_libraries < <(count_read_adt_csv "${project_libraries}" "${adt_library_csv}")

                ADT_index_folder=${project_outs}/$library/adt_index
                ADT_outs=${project_outs}/$library/ADT/

                corrected_fastq=$(realpath -m "${fastq_dirs[0]}/../../KITE_corrected")
                ADT_file="${project_scripts}/adt_files/${ADT_file}.csv"

                if [[ ! -f "$ADT_file" ]]; then
                    echo -e "\033[0;31mERROR: File does not exist: $ADT_file. Check metadata and adt_files folder\033[0m"
                    exit 1
                fi
                
                if [ "${count_submitted}" = "yes" ]; then
                    sbatch_dependency="--dependency=afterok:${job_id}"
                else
                    sbatch_dependency=""
                fi

                sbatch $sbatch_dependency <<EOF
#!/bin/bash
#SBATCH --job-name ${library}_ADT
#SBATCH --output ${project_outs}/logs/kite_${library}.out
#SBATCH --error ${project_outs}/logs/kite_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=96GB
#SBATCH --time=72:00:00

source "${oscar_dir}/functions.sh"

# Log input variables
log "Input variables:"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "count_container: ${count_container}"
log "ADT_index_folder: ${ADT_index_folder}"
log "corrected_fastq: ${corrected_fastq}"
log "ADT_outs: ${ADT_outs}"
log "ADT_file: ${ADT_file}"
log "fastq_dirs: ${fastq_dirs}"
log "fastq_libraries: ${fastq_libraries}"
log "TMPDIR: ${TMPDIR}"
log "sbatch_dependency: ${sbatch_dependency}"

log ""

cd ${project_outs}/${library}

# Create required directories
log "Creating directories"
mkdir -p ${ADT_index_folder}/temp
mkdir -p $corrected_fastq
mkdir -p ${ADT_outs}
check_status "Directory creation"

log ""

# Run featuremap
log "Running featuremap"
apptainer run -B /data ${count_container} featuremap ${ADT_file} \
    --t2g ${ADT_index_folder}/FeaturesMismatch.t2g \
    --fa ${ADT_index_folder}/FeaturesMismatch.fa \
    --header --quiet
check_status "featuremap"

log ""

# Running kallisto index
log "Running kallisto index"
apptainer run -B /data ${count_container} kallisto index \
    -i ${ADT_index_folder}/FeaturesMismatch.idx \
    -k 15 ${ADT_index_folder}/FeaturesMismatch.fa
check_status "kallisto index"

log ""

library_out_name=$(echo "$library" | sed 's/_ATAC/_ADT/')

# Check if files already exist
if [ ! -f "${corrected_fastq}/\${library_out_name}_R1.fastq.gz" ] || [ ! -f "${corrected_fastq}/\${library_out_name}_R2.fastq.gz" ]; then
    log "Running ASAP to KITE conversion"
    apptainer run -B /data ${count_container} asap_to_kite \
        -ff "$fastq_dirs" \
        -sp "$fastq_libraries" \
        -of "${corrected_fastq}/\${library_out_name}" \
        -c $(nproc)
    check_status "ASAP to KITE conversion"
else
    log "KITE converted files already exist, skipping conversion"
fi

log ""

# Running kallisto bus
log "Running kallisto bus"
apptainer run -B /data ${count_container} kallisto bus \
    -i ${ADT_index_folder}/FeaturesMismatch.idx \
    -o ${ADT_index_folder}/temp \
    -x 0,0,16:0,16,26:1,0,0 \
    -t \$(nproc) \
    ${corrected_fastq}/\${library_out_name}*
check_status "kallisto bus"

log ""

log "Extracting ATAC barcodes to ${TMPDIR}/OSCAR/737K-arc-v1.txt"
apptainer exec ${count_container} gunzip \
-c /opt/cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz > \
${TMPDIR}/OSCAR/737K-arc-v1.txt
ATAC_whitelist=${TMPDIR}/OSCAR/737K-arc-v1.txt
check_status "gunzip barcodes"

log ""

log "Running bustools correct"
apptainer run -B /data ${count_container} bustools correct \
    -w \${ATAC_whitelist} \
    ${ADT_index_folder}/temp/output.bus \
    -o ${ADT_index_folder}/temp/output_corrected.bus
check_status "bustools correct"

log ""

# Running bustools sort
log "Running bustools sort"
apptainer run -B /data ${count_container} bustools sort \
    -t \$(nproc) \
    -o ${ADT_index_folder}/temp/output_sorted.bus \
    ${ADT_index_folder}/temp/output_corrected.bus
check_status "bustools sort"

log ""

# Running bustools count
log "Running bustools count"
apptainer run -B /data ${count_container} bustools count \
    -o ${ADT_outs} \
    --genecounts \
    -g ${ADT_index_folder}/FeaturesMismatch.t2g \
    -e ${ADT_index_folder}/temp/matrix.ec \
    -t ${ADT_index_folder}/temp/transcripts.txt \
    ${ADT_index_folder}/temp/output_sorted.bus
check_status "bustools count"

log ""

# Cleanup
log "Cleaning up temporary files"
rm -r ${ADT_index_folder}
check_status "Cleanup"

log ""

log "All processing completed successfully"
EOF
            fi
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX*' "${project_libraries}/${library}.csv"; then
    echo "Processing ${library} as an CITE/GEX run"
        echo ""
        echo "For library $library"
        echo ""
        cat ${project_libraries}/${library}.csv
        echo ""
        echo "cellranger multi --id ${library} --csv ${project_libraries}/${library}.csv --localcores 32"

        # Ask the user if they want to submit the indices for FASTQ generation
        read -p "Is this alright? (Y/N): " choice
        while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -p "Is this alright? (Y/N): " choice
        done
        # Process choices
        if [ "${choice}" = "Y" ] || [ "${choice}" = "y" ]; then
            mkdir -p ${project_outs}/logs/
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output ${project_outs}/logs/cellranger_${library}.out
#SBATCH --error ${project_outs}/logs/cellranger_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=96GB
#SBATCH --time=96:00:00

source "${oscar_dir}/functions.sh"

# Log input variables
log "Input variables:"
log "library: ${library}"
log "project_outs: ${project_outs}"
log "oscar_dir: ${oscar_dir}"
log "count_container: ${count_container}"
log "project_libraries: ${project_libraries}"

log ""

# Run the CellRanger multi command
log "Running cellranger multi"
apptainer run -B /data "${count_container}" cellranger multi \
    --id "${library}" \
    --csv "${project_libraries}/${library}.csv" \
    --localcores "\$(nproc)"
check_status "cellranger multi"

log ""

# Clean up temporary files
log "Cleaning up temporary files"
rm -r ${project_outs}/${library}/SC_MULTI_CS ${project_outs}/${library}/_*
check_status "Cleanup"

log "All processing completed successfully!"
EOF
        fi
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine whether ${library} is a GEX or ATAC run. is 'GEX' or 'ATAC present in its library csv file?"
        exit 1
    fi
        # Reset variables
        job_id=""
        fastq_names=""
        fastq_dirs=""
        extra_arguments=""
        count_submitted=""
        sbatch_dependency=""
        ADT_file=""
        adt_library_csv=""
        ADT_index_folder=""
        corrected_fastq=""
        fastq_libraries=""
        ADT_outs=""
        library_out_name=""

    fi
done
