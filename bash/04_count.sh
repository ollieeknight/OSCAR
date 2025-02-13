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

IFS=',' read -r -a project_ids <<< "${project_id}"

output_project_id="${project_ids[0]}"
output_project_dir="${dir_prefix}/${output_project_id}"
output_project_scripts="${output_project_dir}/${output_project_id}_scripts/"
output_project_libraries=${output_project_scripts}/libraries
output_project_outs="${output_project_dir}/${output_project_id}_outs/"

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

# Pull necessary OSCAR containers
check_and_pull_oscar_containers
count_container=${TMPDIR}/OSCAR/oscar-count_latest.sif

# Take the csv files into a list and remove the .csv suffix
libraries=($(ls "${output_project_libraries}" | awk -F/ '{print $NF}' | awk -F. '{print $1}'))
echo "Libraries:"
for library in "${libraries[@]}"; do
    echo "${library}"
done
exit 1
mkdir -p ${output_project_outs}/

# Iterate over each library file to submit counting jobs
for library in "${libraries[@]}"; do
    # Skip processing lines with 'ADT' in the library name
    if grep -q '.*\(HTO\|ADT\).*' "${output_project_libraries}/${library}.csv"; then
#        echo "Processing ${library} as an ADT/HTO library"
        continue
    elif grep -q '.*\(ATAC\).*' "${output_project_libraries}/${library}.csv"; then
        fastq_names=""
        fastq_dirs=""

        read fastq_names fastq_dirs < <(count_read_csv "${output_project_libraries}" "$library")

        extra_arguments=$(count_check_dogma "${output_project_libraries}" "$library")

        read -p "${library} as an ATAC library; process with cellranger-atac? (Y/N): " choice
        while [[ ! $choice =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -p "Submit ${library} for cellranger-atac count? (Y/N): " choice
        done
        # Process choices
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
            mkdir -p ${output_project_outs}/logs
            # Submit the job to slurm for counting
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output ${output_project_outs}/logs/cellranger_${library}.out
#SBATCH --error ${output_project_outs}/logs/cellranger_${library}.out
#SBATCH --ntasks=64
#SBATCH --mem=128GB
#SBATCH --time=96:00:00

source "${oscar_dir}/functions.sh"

log "OSCAR step 4: cellranger-atac mapping"
log "See https://github.com/ollieeknight/OSCAR for more information"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Cellranger flavour      | cellranger-atac"
log "Library                 | $library"
log "Reference genome        | $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/"
log "FASTQ directory         | $fastq_dirs"
log "FASTQ samples           | $fastq_names"
log "Cores                   | \$(nproc)"
log "Extra arguments         | $extra_arguments"
log "----------------------------------------"

echo ""

cd ${output_project_outs}

# Run cellranger-atac count
log "Running cellranger-atac count"

echo ""
apptainer run -B /data ${count_container} cellranger-atac count \
    --id $library \
    --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ \
    --fastqs $fastq_dirs \
    --sample $fastq_names \
    --localcores \$(nproc) \
    $extra_arguments

echo ""

rm -r ${output_project_outs}/$library/_* ${output_project_outs}/$library/SC_ATAC_COUNTER_CS

log "All processing completed successfully!"

EOF
            )
            count_submitted='YES'
            job_id=$(echo "$job_id" | awk '{print $4}')
        elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
            count_submitted='NO'
        else
            echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
        fi

        if grep -q '.*\(ASAP\).*' "${output_project_libraries}/${library}.csv"; then

           temp_library="${library/_ATAC/}"
            ADT_file=""
            for project_id in "${project_ids[@]}"; do
                project_dir="${dir_prefix}/${project_id}"
                project_scripts="${project_dir}/${project_id}_scripts"
                metadata_file="${project_scripts}/metadata/${metadata_file_name}"
                ADT_file="$(extract_adt_file "$metadata_file" "$temp_library")"
                ADT_file="${project_scripts}/adt_files/${ADT_file}.csv"

                # Check if the ADT file exists
                if [[ -f "${ADT_file}" ]]; then
                    break
                fi
            done

            # Check if the ADT file exists
            if [[ ! -f "${ADT_file}" ]]; then
                echo "ERROR: ${ADT_file} not found."
                exit 1
            fi
            
#            echo "DEBUG: ${ADT_file} found."

            # Determine the correct ADT CSV file name by replacing _ATAC with _ADT
            adt_library_csv="${library/_ATAC/_ADT}.csv"
            adt_library_csv="${output_project_libraries}/${adt_library_csv}"

            # Check if the ADT CSV file exists
            if [[ ! -f "${adt_library_csv}" ]]; then
                echo "ERROR: ${adt_library_csv} not found."
                exit 1
            fi

            read -p "Perform ADT counting? (Y/N): " choice

            while [[ ! $choice =~ ^[YyNn]$ ]]; do
                echo "Invalid input. Please enter Y or N."
                read -p "Perform ADT counting? (Y/N): " choice
            done

            if [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                continue
            elif [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

                read fastq_dirs fastq_libraries < <(count_read_adt_csv "${output_project_libraries}" "${adt_library_csv}")

                ADT_index_folder=${output_project_outs}/$library/adt_index
                ADT_outs=${output_project_outs}/$library/ADT/

                corrected_fastq=$(realpath -m "${fastq_dirs[0]}/../../KITE_corrected")
                
                if [ "${count_submitted}" = "YES" ]; then
                    sbatch_dependency="--dependency=afterok:${job_id}"
                else
                    sbatch_dependency=""
                fi

                sbatch $sbatch_dependency <<EOF
#!/bin/bash
#SBATCH --job-name ${library}_ADT
#SBATCH --output ${output_project_outs}/logs/kite_${library}.out
#SBATCH --error ${output_project_outs}/logs/kite_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=96GB
#SBATCH --time=72:00:00

source "${oscar_dir}/functions.sh"

log "OSCAR step 4: ASAP-seq ADT/HTO mapping"
log "See https://github.com/ollieeknight/OSCAR for more information"

echo ""

library_out_name=$(echo "$library" | sed 's/_ATAC/_ADT/')

apptainer exec ${count_container} gunzip -c /opt/cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz > ${TMPDIR}/OSCAR/737K-cratac-v1.txt.gz

ATAC_whitelist=${TMPDIR}/OSCAR/737K-cratac-v1.txt.gz

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "ADT file                | ${ADT_file}"
log "ADT index folder        | ${ADT_index_folder}"
log "Input FASTQ files       | $fastq_dirs"
log "FASTQ to convert        | $fastq_libraries"
log "Corrected FASTQ name    | \${library_out_name}"
log "Corrected FASTQ output  | ${corrected_fastq}/\${library_out_name}/"
log "Cores                   | \$(nproc)"
log "Barcode whitelist       | \${ATAC_whitelist}"
log "----------------------------------------"

echo ""

cd ${output_project_outs}/${library}

# Create required directories
mkdir -p ${ADT_index_folder}/temp
mkdir -p $corrected_fastq
mkdir -p ${ADT_outs}

# Run featuremap
log "Running featuremap"
apptainer run -B /data ${count_container} featuremap ${ADT_file} \
    --t2g ${ADT_index_folder}/FeaturesMismatch.t2g \
    --fa ${ADT_index_folder}/FeaturesMismatch.fa \
    --header --quiet

echo "" 

# Running kallisto index
log "Running kallisto index"
apptainer run -B /data ${count_container} kallisto index \
    -i ${ADT_index_folder}/FeaturesMismatch.idx \
    -k 15 ${ADT_index_folder}/FeaturesMismatch.fa

echo ""

# Check if files already exist
if [ ! -f "${corrected_fastq}/\${library_out_name}/\${library_out_name}_R1.fastq.gz" ] || [ ! -f "${corrected_fastq}/\${library_out_name}/\${library_out_name}_R2.fastq.gz" ]; then
    log "Running ASAP to KITE conversion"
    apptainer run -B /data ${count_container} asap_to_kite \
        -ff "$fastq_dirs" \
        -sp "$fastq_libraries" \
        -of "${corrected_fastq}/\${library_out_name}" \
        -on \${library_out_name} \
        -c \$(nproc)
else
    log "KITE converted files already exist, skipping conversion"
fi

echo ""

# Running kallisto bus
log "Running kallisto bus"
apptainer run -B /data ${count_container} kallisto bus \
    -i ${ADT_index_folder}/FeaturesMismatch.idx \
    -o ${ADT_index_folder}/temp \
    -x 0,0,16:0,16,26:1,0,0 \
    -t \$(nproc) \
    ${corrected_fastq}/\${library_out_name}/*

echo ""

log "Running bustools correct"
apptainer run -B /data ${count_container} bustools correct \
    -w \${ATAC_whitelist} \
    ${ADT_index_folder}/temp/output.bus \
    -o ${ADT_index_folder}/temp/output_corrected.bus

echo ""

# Running bustools sort
log "Running bustools sort"
apptainer run -B /data ${count_container} bustools sort \
    -t \$(nproc) \
    -o ${ADT_index_folder}/temp/output_sorted.bus \
    ${ADT_index_folder}/temp/output_corrected.bus

echo ""

# Running bustools count
log "Running bustools count"
apptainer run -B /data ${count_container} bustools count \
    -o ${ADT_outs} \
    --genecounts \
    -g ${ADT_index_folder}/FeaturesMismatch.t2g \
    -e ${ADT_index_folder}/temp/matrix.ec \
    -t ${ADT_index_folder}/temp/transcripts.txt \
    ${ADT_index_folder}/temp/output_sorted.bus

echo ""

# Cleanup
rm -r ${ADT_index_folder}

log "All processing completed successfully!"

EOF
            fi
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX*' "${output_project_libraries}/${library}.csv"; then
=        echo ""
        echo "For library $library"
        echo ""
        cat ${output_project_libraries}/${library}.csv
        echo ""

        # Ask the user if they want to submit the indices for FASTQ generation
        read -p "${library} is a GEX or CITE library, process with cellranger multi? (Y/N): " choice
        while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -p "Is this alright? (Y/N): " choice
        done
        # Process choices
        if [ "${choice}" = "Y" ] || [ "${choice}" = "y" ]; then
            mkdir -p ${output_project_outs}/logs/
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output ${output_project_outs}/logs/cellranger_${library}.out
#SBATCH --error ${output_project_outs}/logs/cellranger_${library}.out
#SBATCH --ntasks=64
#SBATCH --mem=128GB
#SBATCH --time=96:00:00

source "${oscar_dir}/functions.sh"

log "OSCAR step 4: cellranger multi mapping"
log "See https://github.com/ollieeknight/OSCAR for more information"

echo ""

log "Input variables:"
log "----------------------------------------"
log "Variable                | Value"
log "----------------------------------------"
log "Library                 | ${library}"
log "Cores                   | $(nproc)"
log "----------------------------------------"

echo ""

# Run the CellRanger multi command
log "Running cellranger multi"
apptainer run -B /data "${count_container}" cellranger multi \
    --id "${library}" \
    --csv "${output_project_libraries}/${library}.csv" \
    --localcores "\$(nproc)"

echo ""

# Clean up temporary files
rm -r ${output_project_outs}/${library}/SC_MULTI_CS ${output_project_outs}/${library}/_*

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
