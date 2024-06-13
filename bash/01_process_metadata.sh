#/bin/bash

# Define default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
dir_prefix="${HOME}/scratch/ngs"

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-id)
      project_id="$2"
      shift 2
      ;;
    --dir_prefix)
      dir_prefix="$2"
      shift 2
      ;;
    *)      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

# Check if project_id is empty
if [ -z "$project_id" ]; then
    echo -e "\033[0;31mERROR:\033[0m Please provide a project_id using the --project-id option."
    exit 1
fi

read_length=$(awk -F '"' '/<Read Number="1"/ {print $4}' ${dir_prefix}/${project_id}/${project_id}_bcl/RunInfo.xml)
if [ "${read_length}" -gt 45 ]; then
    run_type="ATAC"
elif [ "${read_length}" -lt 45 ]; then
    run_type="GEX"
else
    echo -e "\033[0;31mERROR:\033[0m Cannot determine run type, please check ${project_dir}/${project_id}_bcl/RunInfo.xml"
    exit 1
fi

echo ""
echo -e "\033[34mINFO:\033[0m ${project_id} is an ${run_type} run with an R1 read length of ${read_length}"

# Define project directory using the dir_prefix
project_dir="${dir_prefix}/${project_id}"
project_scripts="${project_dir}/${project_id}_scripts"

# In case indices folder is present, remove indices folder to start fresh
if [ -d "${project_scripts}/indices" ]; then
    rm -r "${project_scripts}/indices"
fi

mkdir -p ${project_scripts}/indices
project_indices=${project_scripts}/indices

metadata_file=${project_scripts}/metadata/metadata.csv

if [[ ! -f "${metadata_file}" ]]; then
    echo -e "\033[0;31mERROR:\033[0m Metadata file not found for ${project_id}"
    exit 1
fi

while IFS=',' read -r assay experiment_id historical_number replicate modality chemistry index_type index species n_donors adt_file; do
    # Skip the first header line
    if [[ ${assay} == assay ]]; then
        continue
    fi

    # Add some lines for output readability
    echo ""
    echo "-------------"
    echo ""
    # Print the line being processed
    echo "Processing metadata line: ${line}"

    if [ "${chemistry}" != "NA" ]; then
        # Create the output file with ${chemistry} included
        sample="${project_indices}/${assay}_${index_type}_${modality}_${chemistry}"
	echo $sample
        output_file="${sample}.csv"
    else
        # Create the output file without $chemistry
        sample="${project_indices}/${assay}_${index_type}_${modality}"
	echo $sample
        output_file="${sample}.csv"
    fi

    # Check if the csv file already exists
    if [ ! -f "${output_file}" ]; then
        # If the file doesn't exist, create it and add the header and sample
        echo "Output file ${output_file} does not exist, creating csv and appending ${assay}_${experimental_id}_exp${historical_number}_lib${replicate}_${modality}"
        echo "lane,sample,index" > "${output_file}"
        echo "*,${assay}_${experimental_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "${output_file}"
    else
        # If the file exists, just add sample
        echo "Output file ${output_file} already exists, appending appending ${assay}_${experimental_id}_exp${historical_number}_lib${replicate}_${modality}"
        echo "*,${assay}_${experimental_id}_exp${historical_number}_lib${replicate}_${modality},${index}" >> "$output_file"
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
    bash ${oscar_dir}/02_fastq.sh --project-id ${project_id}
else
    continue
fi
