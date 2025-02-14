#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")  # Get the directory of the current script
source "${oscar_dir}/functions.sh"  # Source the functions script
dir_prefix="${HOME}/scratch/ngs"  # Default directory prefix
metadata_file_name="metadata.csv"  # Default metadata file name

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

check_project_id  # Function to check the project ID

# Set project directories
project_dir="${dir_prefix}/${project_id}"
project_scripts="${project_dir}/${project_id}_scripts"
project_indices="${project_scripts}/indices"

# Set metadata file path
metadata_file="${project_scripts}/metadata/${metadata_file_name}"
check_metadata_file "${metadata_file}"  # Function to check the metadata file

# Determine the run type
run_type=$(check_run_type "${project_id}" "${dir_prefix}")

echo ""
echo -e "\033[34mINFO:\033[0m Detected a(n) ${run_type} run"

# Recreate indices folder
if [ -d "${project_indices}" ]; then
  rm -r "${project_indices}"  # Remove existing indices folder
fi

mkdir -p "${project_indices}"  # Create new indices folder

# Read and process the metadata file
while IFS=',' read -r assay experiment_id historical_number replicate modality chemistry index_type index species n_donors adt_file || [[ -n "$assay" ]]; do
  # Skip the first header line
  if [[ "${assay}" != "assay" ]]; then
    echo ""
    echo "Processing library component ${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality}"

    # Determine the output file name based on chemistry
    if [ "${chemistry}" != "NA" ] && ( [ "${assay}" == "CITE" ] || [ "${assay}" == "GEX" ] ); then
      sample="${project_indices}/${assay}_${index_type}_${modality}_${chemistry}"
    else
      sample="${project_indices}/${assay}_${index_type}_${modality}"
    fi

    output_file="${sample}.csv"

    # Check if the output file exists
    if [ ! -f "${output_file}" ]; then
      echo "Creating ${output_file} and appending line"
      echo "lane,sample,index" > "${output_file}"  # Create new CSV file with header
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"  # Add data line
    else
      echo "Appending line to ${output_file}"
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"  # Append data line
    fi
  fi
done < "${metadata_file}"

# Ask the user if they want to submit the indices for FASTQ generation
echo ""
echo -e "\033[0;33mINPUT:\033[0m Would you like to proceed to FASTQ demultiplexing? (y/n)"
echo -e "bash ${oscar_dir}/02_fastq.sh --project-id ${project_id}"
read -r choice
while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
  echo "Invalid input. Please enter y or n"
  read -r choice
done

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
  bash "${oscar_dir}/02_fastq.sh" --project-id "${project_id}"  # Run the FASTQ script
else
  exit 0  # Exit the script
fi
