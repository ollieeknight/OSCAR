#! /bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"
gene_expression_options=""
vdj_options=""
adt_options=""
mode=""

# Function to display help message
display_help() {
        echo "Usage: $0 [options]"
        echo ""
        echo "Options:"
        echo "  --project-id <id>                             Set the project ID. Can be multiple (comma-separated)"
        echo "  --dir-prefix <path>                           Set the directory prefix (default: ${HOME}/scratch/ngs)"
        echo "  --gene-expression-options <id1,id2>           Define options for gene expression processing (comma-separated)"
        echo "  --vdj-options <id1,id2>                       Define options for VDJ processing (comma-separated)"
        echo "  --adt-options <id1,id2>                       Define options for ADT processing (comma-separated)"
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

# Split project_ids, gene_expression_options, vdj_options, and adt_options into arrays
IFS=',' read -r -a project_ids <<< "${project_id}"
IFS=';' read -r -a gene_expression_options <<< "${gene_expression_options}"
IFS=';' read -r -a vdj_options <<< "${vdj_options}"
IFS=';' read -r -a adt_options <<< "${adt_options}"

# Define output directories
output_project_id="${project_ids[0]}"
output_project_dir="${dir_prefix}/${output_project_id}"
output_project_scripts="${output_project_dir}/${output_project_id}_scripts/"
output_project_libraries=${output_project_scripts}/libraries

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

# List index files and extract flowcell ID from RunInfo.xml for the first project ID
flowcell_id=$(grep "<Flowcell>" "${project_dir}/${project_ids[0]}_bcl/RunInfo.xml" | sed -e 's|.*<Flowcell>\(.*\)</Flowcell>.*|\1|')

# Pull necessary OSCAR containers
check_and_pull_oscar_containers

# Validate the mode
validate_mode "${mode}"

# Print options for gene expression, VDJ, and ADT if mode is GEX
print_options "${mode}" "Gene expression" "${gene_expression_options[@]}"
print_options "${mode}" "VDJ-B/T" "${vdj_options[@]}"
print_options "${mode}" "ADT/HTO" "${adt_options[@]}"

# Check if the libraries folder already exists, and remove it if it does
if [ -d "${output_project_libraries}" ]; then
        rm -r "${output_project_libraries}"
fi

# Create the libraries folder
mkdir -p "${output_project_libraries}"

# Process each project ID
for project_id in "${project_ids[@]}"; do
        project_dir=${dir_prefix}/${project_id}
        project_scripts=${project_dir}/${project_id}_scripts

        # Define the metadata file path based on the project_id
        metadata_file="${project_scripts}/metadata/metadata.csv"

        # Check that the metadata file is available
        if [ ! -f "${metadata_file}" ]; then
                echo -e "\033[0;31mERROR:\033[0m Metadata file for ${project_id} not found, please check path"
                exit 1
        fi

        # Determine the run type
        run_type=$(check_run_type "${project_id}" "${dir_prefix}")

        echo -e "\033[34mINFO:\033[0m ${project_id} is an ${run_type} run, processing appropriately"

        # Iterate through each line in metadata.csv
        while IFS=',' read -r assay experiment_id historical_number replicate modality chemistry index_type index species n_donors adt_file; do
                # Skip the header line
                if [[ ${assay} == "assay" ]]; then
                        continue
                fi

                # Define the library name
                library="${assay}_${experiment_id}_exp${historical_number}_lib${replicate}"

                # Determine the full modality
                full_modality=$(determine_full_modality "${modality}" "${library}")
                if [ ${?} -eq 1 ]; then
                        continue
                fi

                # Define the output file for the library
                library_output=${output_project_libraries}/${library}.csv

                # Handle different run types and modes
                if [[ ${run_type} == 'GEX' && ${mode} == "GEX" ]]; then
                        handle_gex_mode
                elif [[ ${run_type} == 'ATAC' && ${mode} == "GEX" ]]; then
                        handle_atac_mode
                elif [[ ${run_type} == 'ATAC' && ${mode} == "ATAC" ]]; then
                        handle_atac_mode
                fi

        done < "${metadata_file}"
done

# Ask the user if they want to submit the libraries for counting
echo "Would you like to proceed to counting? (Y/N)"
echo -e "bash ${oscar_dir}/04_count.sh --project-id ${output_project_id}"
read -r choice
while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
        echo "Invalid input. Please enter y or n"
        read -r choice
done

# Process choices
if [ "${choice}" = "Y" ] || [ "${choice}" = "y" ]; then
        bash ${oscar_dir}/04_count.sh --project-id ${output_project_id}
elif [ "${choice}" = "N" ] || [ "${choice}" = "n" ]; then
        :
else
        echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
fi
