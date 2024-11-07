#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")  # Get the directory of the current script
source "${oscar_dir}/functions.sh"  # Source the functions script
dir_prefix="${HOME}/scratch/ngs"  # Default directory prefix
metadata_file_name="metadata.csv"  # Default metadata file name

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
  if [[ "$1" == --* ]]; then
    var_name=$(echo "$1" | sed 's/--//; s/-/_/')  # Convert --option-name to option_name
    declare "$var_name"="$2"  # Declare the variable with the given value
    shift 2  # Shift to the next pair of arguments
  else
    echo "Invalid option: $1"  # Print error message for invalid option
    exit 1  # Exit with error code
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
echo -e "\033[34mINFO:\033[0m Detected a(n) ${run_type} run"

# Recreate indices folder
if [ -d "${project_indices}" ]; then
  rm -r "${project_indices}"  # Remove existing indices folder
fi

mkdir -p "${project_indices}"  # Create new indices folder

# Read and process the metadata file
while IFS=',' read -r assay experiment_id historical_number replicate modality chemistry index_type index species n_donors adt_file; do
  # Skip the first header line
  if [[ "${assay}" != "assay" ]]; then
    echo ""
    echo "Processing metadata line"
    echo "${assay},${experiment_id},${historical_number},${replicate},${modality},${chemistry},${index_type},${index},${species},${n_donors},${adt_file}"

    # Determine the output file name based on chemistry
    if [ "${chemistry}" != "NA" ]; then
      sample="${project_indices}/${assay}_${index_type}_${modality}_${chemistry}"
    else
      sample="${project_indices}/${assay}_${index_type}_${modality}"
    fi

    output_file="${sample}.csv"

    # Check if the output file exists
    if [ ! -f "${output_file}" ]; then
      echo "Output file ${output_file} does not exist, creating csv"
      echo "lane,sample,index" > "${output_file}"  # Create new CSV file with header
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"  # Add data line
    else
      echo "Output file ${output_file} already exists, appending"
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"  # Append data line
    fi
  fi
done < "${metadata_file}"

# Ask the user if they want to submit the indices for FASTQ generation
echo ""
echo -e "\033[0;33mINPUT REQUIRED:\033[0m Would you like to proceed to FASTQ demultiplexing? (y/n)"
read -r choice
while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter y or n"
    read -r choice
done

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
    echo "Submitting: bash ${oscar_dir}/02_fastq.sh --project-id ${project_id}"
    bash "${oscar_dir}/02_fastq.sh" --project-id "${project_id}"  # Run the FASTQ script
else
    exit 0  # Exit the script
fi
