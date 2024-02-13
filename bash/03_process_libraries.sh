#!/bin/bash

# usage is: bash 03_process_libraries.sh --project-id 'project_id' --gene-expression-options option1,setting1;option2,setting 2 --vdj-options option1,setting1;option2,setting 2 --adt-options option1,setting1;option2,setting 2
# such as bash 03_process_libraries.sh --project-id K002,K001 --prefix $HOME/scratch/ngs/
# The output project id should be the first! Above, output libraries will go into K002.

# make sure your processing folder structure is similar to follows:

# |── {project_id}_bcl/
# │   ├── Config/
# │   ├── CopyComplete.txt
# │   ├── Data/
# │   ├── InterOp/
# │   ├── Logs/
# │   ├── Recipe/
# │   ├── RTA3.cfg
# │   ├── RTAComplete.txt
# │   ├── RunInfo.xml
# │   ├── RunParameters.xml
# │   ├── SequenceComplete.txt
# │   └── Thumbnail_Images/
# |── {project_id}_fastq/
# │   ├── FASTQ_1/
# │   ├── FASTQ_2/
# │   ├── ...
# │   ├── FASTQ_n/
# └── {project_id}_scripts/
#     ├── adt_files/
#     ├── indices/
#     └── metadata/

# Define default values
OSCAR_script_dir=$(dirname "${BASH_SOURCE[0]}")
prefix="$HOME/scratch/ngs"

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-id)
      project_ids_in="$2"
      shift 2
      ;;
    --prefix)
      prefix="$2"
      shift 2
      ;;
    --gene-expression-options)
      gene_expression_options="$2"
      shift 2
      ;;
    --vdj-options)
      vdj_options="$2"
      shift 2
      ;;
    --adt-options)
      adt_options="$2"
      shift 2
      ;;
    *)      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

IFS=',' read -r -a project_ids <<< "$project_ids_in"
IFS=';' read -r -a gene_expression_options <<< "$gene_expression_options"
IFS=';' read -r -a vdj_options <<< "$vdj_options"
IFS=';' read -r -a adt_options <<< "$adt_options"

main_project_id="${project_ids[0]}"

# Check if project_id is empty
if [ -z "${project_ids[0]}" ]; then
  echo ""
  echo "Please provide at least one project_id using the --project-id option"
  echo ""
  echo "Option fields can be left blank, and you can find options here https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts"
  echo ""
  exit 1
fi

# Check if gene_expression_options array has entries
if [ ${#gene_expression_options[@]} -gt 0 ]; then
    echo "Gene expression options set as:"
    for option in "${gene_expression_options[@]}"; do
        echo "$option"
    done
else
    echo "No options set for gene expression"
fi

# Check if vdj_options array has entries
if [ ${#vdj_options[@]} -gt 0 ]; then
    echo "VDJ-B/T options set as:"
    for option in "${vdj_options[@]}"; do
        echo "$option"
    done
else
    echo "No options set for VDJ-B/T"
fi

# Check if adt_options array has entries
if [ ${#adt_options[@]} -gt 0 ]; then
    echo "ADT/HTO options set as:"
    for option in "${adt_options[@]}"; do
        echo "$option"
    done
else
    echo "No options set for ADT/HTO"
fi

output_folder="${prefix}/${main_project_id}/${main_project_id}_scripts/libraries"

# Check if the libraries folder already exists, and remove it if it does
if [ -d "$output_folder" ]; then
  rm -r "$output_folder"
fi

mkdir -p "$output_folder"

for project_id in "${project_ids[@]}"; do
    echo ""
    echo "Processing project_id: $project_id"
    echo ""

	project_dir=$prefix/$project_id
    script_dir=${project_dir}/${project_id}_scripts

    # Define the metadata file path based on the project_id
    metadata_file="$project_dir/${project_id}_scripts/metadata/metadata.csv"

    # Check that the singularity container is available
    if [ ! -f "${metadata_file}" ]; then
        echo "Metadata file for $project_id not found, please check path"
        exit 1
    fi

    # Iterate through each line in metadata.csv
    while IFS= read -r line; do
        # Skip the header line
        if [[ $line == assay* ]]; then
            continue
        fi

        # Split the line into fields
        IFS=',' read -r -a fields <<< "$line"

        # Assign field values to variables
        assay="${fields[0]}"
        experiment_id="${fields[1]}"
        historical_number="${fields[2]}"
        replicate="${fields[3]}"
        modality="${fields[4]}"
        chemistry="${fields[5]}"
        species="${fields[8]}"
        adt_file="${fields[10]}"
        library="${assay}_${experiment_id}_exp${historical_number}_lib${replicate}"

        # Define the library_output path
        library_output=${output_folder}/${library}.csv
        # Check if the assay is not ASAP or ATAC; essentially is it GEX, CITE, or a MULTIOME or DOGMA gene expression run
        if [[ $assay != 'ASAP' && $modality != 'ATAC' ]]; then
            # Check if the sample library already exists
            if [ ! -f "$library_output" ]; then
            # What species is specified in the metadata file?
                # Is it a human run?
                if [[ "$species" =~ ^(Human|human|Hs|hs)$ ]]; then
                    echo "Writing human reference files for $library"
                    echo "[gene-expression]" >> "${library_output}"
                    echo "reference,/fast/work/groups/ag_romagnani/ref/hs/GRCh38-hardmasked-optimised-arc" >> "${library_output}"
                        # Add options if there are gene expression-specific options specified by --gene-expression-options
                        if [ -n "$gene_expression_options" ] && [ "$gene_expression_options" != "NA" ]; then
                                IFS=';' read -ra values <<< "$gene_expression_options"
                                for value in "${values[@]}"; do
                                        echo "$value" >> "${library_output}"
                                done
                        fi
                    # If this is a DOGMA-seq or MULTIOME run, to specify chemistry
                    if [ "$assay" == "DOGMA" ] || [ "$assay" == "MULTIOME" ]; then
                        echo "chemistry,ARC-v1" >> "${library_output}"
                    fi
                    echo "" >> "${library_output}"
                    echo "[vdj]" >> "${library_output}"
                    echo "reference,/fast/work/groups/ag_romagnani/ref/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0" >> "${library_output}"
                    # Add options if there are VDJ-specific options specified by --vdj-options
                    if [ "$vdj_options" != "NA" ]; then
                        IFS=',' read -ra values <<< "$vdj_options"
                        for value in "${values[@]}"; do
                            echo "$value" >> "${library_output}"
                        done
                    fi
                # Is it a mouse run?
                elif [[ "$species" =~ ^(Mouse|mouse|Mm|mm)$ ]]; then
                    echo "[gene-expression]" >> "${library_output}"
                    echo "reference,/fast/work/groups/ag_romagnani/ref/mm/mouse_mm10_optimized_reference_v2" >> "${library_output}"
                    # Add options if there are gene expression-specific options specified by --gene-expression-options
                        if [ -n "$gene_expression_options" ] && [ "$gene_expression_options" != "NA" ]; then
                                IFS=';' read -ra values <<< "$gene_expression_options"
                                for value in "${values[@]}"; do
                                        echo "$value" >> "${library_output}"
                                done
                        fi
                    # If this is a DOGMA-seq or MULTIOME run, to specify chemistry
                    if [ "$assay" == "DOGMA" ] || [ "$assay" == "MULTIOME" ]; then
                        echo "chemistry,ARC-v1" >> "${library_output}"
                    fi
                    echo "" >> "${library_output}"
                    echo "[vdj]" >> "${library_output}"
                    echo "reference,/fast/work/groups/ag_romagnani/ref/mm/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0" >> "${library_output}"
                    # Add options if there are VDJ-specific options specified by --vdj-options
                    if [ -n "$vdj_options" ] && [ "$vdj_options" != "NA" ]; then
                        IFS=',' read -ra values <<< "$vdj_options"
                        for value in "${values[@]}"; do
                            echo "$value" >> "${library_output}"
                        done
                    fi
                fi
                # Is there paired ADT data for this sample?
                if [ "$adt_file" != "NA" ]; then
                    echo "" >> "${library_output}"
                    echo "[feature]" >> "${library_output}"
                    echo "reference,$script_dir/adt_files/${adt_file}.csv" >> "${library_output}"
                    # Add options if there are ADT/HTO-specific options specified by --adt-options
                    if [ "$adt_options" != "" ]; then
                        IFS=',' read -ra values <<< "$adt_options"
                        for value in "${values[@]}"; do
                            echo "$value" >> "${library_output}"
                        done
                    fi
                else
                    :
                fi
                echo "Writing ${modality} for ${library}"
                echo "" >> "${library_output}"
                echo "[libraries]" >> "${library_output}"
                echo "fastq_id,fastqs,feature_types" >> "${library_output}"
            fi
            # Determine the full-length modality name based on its shortened name
            if [ "$modality" = "GEX" ]; then
                full_modality='Gene Expression'
            elif [ "$modality" = "ADT" ] || [ "$modality" = "HTO" ]; then
                full_modality='Antibody Capture'
            elif [ "$modality" = "VDJ-T" ]; then
               full_modality='VDJ-T'
            elif [ "$modality" = "VDJ-B" ]; then
                full_modality='VDJ-B'
            elif [ "$modality" = "CRISPR" ]; then
                full_modality='CRISPR Guide Capture'
            fi
            # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
            declare -A unique_lines
            # Recursively search for FASTQ files in the project_id FASTQ folder
            for folder in "$project_dir/${project_id}_fastq"/*/outs; do
                matching_fastq_files=($(find "$folder" -type f -name "${library}*${modality}*" | sort -u))
                for fastq_file in "${matching_fastq_files[@]}"; do
                    # Extract the directory containing the FASTQ file
                    directory=$(dirname "$fastq_file")
                    # Extract the modified name from the FASTQ file
                    fastq_name=$(basename "$fastq_file" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                    # Create a unique identifier the FASTQ file
                    line_identifier="$fastq_name,$directory,$full_modality"
                    # Check if the line has already been added
                    if [ ! -v unique_lines["$line_identifier"] ]; then
                        unique_lines["$line_identifier"]=1
                        echo "$fastq_name,$directory,$full_modality" >> "${library_output}"
                        echo "Writing $fastq_name,$directory,$full_modality to $library"
                    fi
                done
            done
        # Check if the assay is an ATAC (ASAP) or DOGMA ATAC run
        elif [[ $modality == 'ATAC' ]] || [[ $assay == 'DOGMA' ]]; then
            # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
            declare -A unique_lines
            # Recursively search for files in the FASTQ folder
            for folder in "$project_dir/${project_id}_fastq"/*/outs; do
                matching_fastq_files=($(find "$folder" -type f -name "${library}*${modality}*" | sort -u))
                for fastq_file in "${matching_fastq_files[@]}"; do
                    # Extract the directory containing the fastq file
                    directory=$(dirname "$fastq_file")
                    # Extract the modified name from the fastq file
                    fastq_name=$(basename "$fastq_file" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                    # Create a unique identifier for FASTQ file
                    line_identifier="{$library}_${modality}"
                    # Check if the line has already been added
                    if [ ! -v unique_lines["$line_identifier"] ]; then
                        unique_lines["$line_identifier"]=1
                        echo "Adding $directory for ${library}_${modality}"
                        echo "$directory" >> "${library_output}"
                        echo "Writing $directory to $library"
                    fi
                done
            done
        # By process of exclusion, only remaining sample lines should be the ADT or HTO from an ASAP-seq run
        elif [[ $modality == 'ADT' ]] || [[ $modality == 'HTO' ]]; then
            echo "This should be an ASAP ADT file and will be processed later"
        fi
    echo ""
    echo "-------------"
    echo ""
    done < "$metadata_file"
done

# Ask the user if they want to submit the indices for FASTQ generation
echo "Would you like to submit the next step to perform read counting? (Y/N)"
read -r choice

# Process choices
if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
    echo "Submitting: bash ${OSCAR_script_dir}/04_count.sh --project-id ${project_id}"
    bash ${OSCAR_script_dir}/04_count.sh --project-id ${project_id}
elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
    :
else
    echo "Invalid choice. Exiting"
fi
