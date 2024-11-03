#!/bin/bash

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

check_project_id

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
  if [[ "$1" == --* ]]; then
    # Remove leading dashes and replace hyphens with underscores
    var_name=$(echo "$1" | sed 's/--//; s/-/_/')
    declare "$var_name"="$2"
    shift 2
  else
    echo "Invalid option: $1"
    exit 1
  fi
done

project_dir="${dir_prefix}/${project_id}"
project_scripts="${project_dir}/${project_id}_scripts"
project_indices="${project_scripts}/indices"

metadata_file="${project_scripts}/metadata/${metadata_file_name}"
check_metadata_file "${metadata_file}"

run_type=$(check_run_type "${project_id}" "${dir_prefix}")
echo -e "\033[34mINFO:\033[0m Detected a(n) ${run_type} run"

# Recreate indices folder
if [ -d "${project_indices}" ]; then
  rm -r "${project_indices}"
fi

mkdir -p "${project_indices}"

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
      echo "lane,sample,index" > "${output_file}"
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"
    else
      echo "Output file ${output_file} already exists, appending"
      echo "*,${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"
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
    bash "${oscar_dir}/02_fastq.sh" --project-id "${project_id}"
else
    exit 0
fi
