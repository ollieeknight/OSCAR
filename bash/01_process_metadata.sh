#!/bin/bash

# usage is: bash 01_process_metadata.sh --project-id project_id --prefix $HOME/path/before/
# such as bash 01_process_metadata.sh --project-id K001 --prefix $HOME/scratch/ngs/

# make sure your processing folder structure is similar to follows:

# |── {project_id}_bcl/
# │   ├── Config/
# │   ├── CopyComplete.txt
# │   ├── Data/
# │   ├── InterOp/
# │   ├── Logs/
# │   ├── Recipe/
# │   ├── RTA3.cfg
# │   ├── RTAComplete.txt
# │   ├── RunInfo.xml
# │   ├── RunParameters.xml
# │   ├── SequenceComplete.txt
# │   └── Thumbnail_Images/
# └── {project_id}_scripts/
#     ├── adt_files/
#     └── metadata/ # METADATA MUST BE IN THIS FOLDER!

# Define default values
OSCAR_script_dir=$(dirname "${BASH_SOURCE[0]}")
prefix="$HOME/scratch/ngs"

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-id)
      project_id="$2"
      shift 2
      ;;
    --prefix)
      prefix="$2"
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

read_length=$(awk -F '"' '/<Read Number="1"/ {print $4}' ${prefix}/${project_id}/${project_id}_bcl/RunInfo.xml)
if [ "$read_length" -gt 45 ]; then
    run_type="ATAC"
elif [ "$read_length" -lt 45 ]; then
    run_type="GEX"
else
    echo -e "\033[0;31mERROR:\033[0m Cannot determine run type, please check ${project_dir}/${project_id}_bcl/RunInfo.xml"
    exit 1
fi

echo ""
echo -e "\033[34mINFO:\033[0m $project_id is an $run_type run, processing appropriately"

# Define project directory using the prefix
project_dir="${prefix}/$project_id"
script_dir="${project_dir}/${project_id}_scripts"

# In case indices folder is present, remove indices folder to start fresh
if [ -d "$script_dir/indices" ]; then
    rm -r "$script_dir/indices"
fi

mkdir -p $script_dir/indices
indices_folder=$script_dir/indices

metadata_file=$script_dir/metadata/metadata.csv

if [[ ! -f "$metadata_file" ]]; then
    echo -e "\033[0;31mERROR:\033[0m Metadata file not found for $project_id"
    exit 1
fi

# Iterate through each sample sub-library in metadata.csv
while IFS= read -r line; do
    # Skip the first header line
    if [[ $line == assay* ]]; then
        continue
    fi
    # Add some lines for output readability
    echo ""
    echo "-------------"
    echo ""
    # Print the line being processed
    echo "Processing metadata line: $line"
    # Split the line into fields
    IFS=',' read -r -a fields <<< "$line"

    # Debug prints
    assay="${fields[0]}"
    experimental_id="${fields[1]}"
    historical_id="${fields[2]}"
    replicate="${fields[3]}"
    modality="${fields[4]}"
    chemistry="${fields[5]}"
    index_type="${fields[6]}"
    index="${fields[7]}"
    if [ "$chemistry" != "NA" ]; then
        # Create the output file with $chemistry included
        sample="$indices_folder/${assay}_${index_type}_${modality}_${chemistry}"
        output_file="${sample}.csv"
    else
        # Create the output file without $chemistry
        sample="$indices_folder/${assay}_${index_type}_${modality}"
        output_file="${sample}.csv"
    fi
    # Check if the csv file already exists
    if [ ! -f "$output_file" ]; then
        # If the file doesn't exist, create it and add the header and sample
        echo "Output file $output_file does not exist, creating csv and appending ${assay}_${experimental_id}_exp${historical_id}_lib${replicate}_${modality}"
        echo "lane,sample,index" > "$output_file"
        echo "*,${assay}_${experimental_id}_exp${historical_id}_lib${replicate}_${modality},${index}" >> "$output_file"
    else
        # If the file exists, just add sample
        echo "Output file $output_file already exists, appending appending ${assay}_${experimental_id}_exp${historical_id}_lib${replicate}_${modality}"
        echo "*,${assay}_${experimental_id}_exp${historical_id}_lib${replicate}_${modality},${index}" >> "$output_file"
    fi
done < "$metadata_file"

# Ask the user if they want to submit the indices for FASTQ generation
echo ""
echo -e "\033[0;33mINPUT REQUIRED:\033[0m Would you like to proceed to FASTQ demultiplexing? (Y/N)"
read -r choice

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
    echo "Submitting: bash ${OSCAR_script_dir}/02_fastq.sh --project-id ${project_id}"
    bash ${OSCAR_script_dir}/02_fastq.sh --project-id ${project_id}
elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
else
    echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
fi
