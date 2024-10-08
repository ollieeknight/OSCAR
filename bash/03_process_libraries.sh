#!/bin/bash

# Define default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
dir_prefix="${HOME}/scratch/ngs"
gene_expression_options=""
vdj_options=""
adt_options=""
mode=""

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-id)
      project_ids="$2"
      shift 2
      ;;
    --dir-prefix)
      dir_prefix="$2"
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
    --mode)
      mode="$2"
      shift 2
      ;;
    *)
      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

IFS=',' read -r -a project_ids <<< "${project_ids}"
IFS=';' read -r -a gene_expression_options <<< "${gene_expression_options}"
IFS=';' read -r -a vdj_options <<< "${vdj_options}"
IFS=';' read -r -a adt_options <<< "${adt_options}"

output_project_id="${project_ids[0]}"

# Check if mode is specified
if [[ -z "${mode}" ]]; then
    echo -e "\033[0;31mERROR:\033[0m  Please specify the mode using the --mode option (ATAC or GEX)."
    exit 1
fi

# Validate mode
if [[ "${mode}" != "ATAC" && "${mode}" != "GEX" ]]; then
    echo -e "\033[0;31mERROR:\033[0m Invalid mode. Mode must be either 'ATAC' or 'GEX'."
    exit 1
fi

# Check if project_id is empty
if [ -z "${project_ids[0]}" ]; then
    echo ""
    echo -e "\033[0;31mERROR:\033[0m Please provide at least one project_id using the --project-id option"
    echo ""
    echo -e "\033[0;31mERROR:\033[0m Option fields can be left blank, and you can find options here https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts"
    echo ""
    exit 1
fi

echo ""

# Check if gene_expression_options array has entries
if [ ${#gene_expression_options[@]} -gt 0 ]; then
    echo "Gene expression options set as:"
    for option in "${gene_expression_options[@]}"; do
        echo "${option}"
    done
else
    echo "No options set for gene expression"
fi

# Check if vdj_options array has entries
if [ ${#vdj_options[@]} -gt 0 ]; then
    echo "VDJ-B/T options set as:"
    for option in "${vdj_options[@]}"; do
        echo "${option}"
    done
else
    echo "No options set for VDJ-B/T"
fi

# Check if adt_options array has entries
if [ ${#adt_options[@]} -gt 0 ]; then
    echo "ADT/HTO options set as:"
    for option in "${adt_options[@]}"; do
        echo "${option}"
    done
else
    echo "No options set for ADT/HTO"
fi

output_project_scripts="${dir_prefix}/${output_project_id}/${output_project_id}_scripts/"

# Check if the libraries folder already exists, and remove it if it does
if [ ! -d "${output_project_scripts}" ]; then
    echo -e "\033[0;31mERROR:\033[0m Scripts folder does not exist, did you enter a wrong project ID?"
    exit 1
fi

output_project_libraries=${output_project_scripts}/libraries

# Check if the libraries folder already exists, and remove it if it does
if [ -d "${output_project_libraries}" ]; then
  rm -r "${output_project_libraries}"
fi

mkdir -p "${output_project_libraries}"

for project_id in "${project_ids[@]}"; do
    echo ""
    project_dir=${dir_prefix}/${project_id}
    project_scripts=${project_dir}/${project_id}_scripts

    # Define the metadata file path based on the project_id
    metadata_file="${project_scripts}/metadata/metadata.csv"

    # Check that the singularity container is available
    if [ ! -f "${metadata_file}" ]; then
        echo -e "\033[0;31mERROR:\033[0m Metadata file for ${project_id} not found, please check path"
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

    echo -e "\033[34mINFO:\033[0m ${project_id} is an ${run_type} run, processing appropriately"

    # Iterate through each line in metadata.csv
    while IFS=',' read -r assay experiment_id historical_number replicate modality chemistry index_type index species n_donors adt_file; do
        # Skip the header line
        if [[ ${assay} == "assay" ]]; then
            continue
        fi

        echo ""
        echo "-------------"
        echo ""

        library="${assay}_${experiment_id}_exp${historical_number}_lib${replicate}"

        # Determine the full-length modality name based on its shortened name
        if [ "${modality}" = "GEX" ]; then
            full_modality='Gene Expression'
        elif [ "${modality}" = "ADT" ] || [ "${modality}" = "HTO" ]; then
            full_modality='Antibody Capture'
        elif [ "${modality}" = "VDJ-T" ]; then
           full_modality='VDJ-T'
        elif [ "${modality}" = "VDJ-B" ]; then
            full_modality='VDJ-B'
        elif [ "${modality}" = "CRISPR" ]; then
            full_modality='CRISPR Guide Capture'
        elif [ "${modality}" = "GENO" ]; then
            continue
        fi

        echo "Adding ${full_modality} for ${library}"

        # Define the library_output path
        library_output=${output_project_libraries}/${library}.csv

        # If this run is of ATAC samples, alongside ADT/HTO
        if [[ ${run_type} == 'GEX' && ${mode} == "GEX" ]]; then
            if [[ "${modality}" == 'GEX' ]]; then
                # Check if the sample library already exists
                if [ ! -f "${library_output}" ]; then
                    # What species is specified in the metadata file?
                    # Is it a human run?
                    if [[ "${species}" =~ ^(Human|human|Hs|hs)$ ]]; then
                        echo "Writing human reference files for ${library}"
                        echo "[gene-expression]" >> "${library_output}"
                        echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-hardmasked-optimised-arc" >> "${library_output}"
#                       echo "probe-set,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/frp-probes/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv" >> "${library_output}"
                        echo "create-bam,true" >> "${library_output}"
                        # Add options if there are gene expression-specific options specified by --gene-expression-options
                        if [ -n "${gene_expression_options}" ] && [ "${gene_expression_options}" != "NA" ]; then
                                IFS=';' read -ra values <<< "${gene_expression_options}"
                                for value in "${values[@]}"; do
                                    echo "${value}" >> "${library_output}"
                                done
                        fi
                        # If this is a DOGMA-seq or MULTIOME run, to specify chemistry
                        if [ "${assay}" == "DOGMA" ] || [ "${assay}" == "MULTIOME" ]; then
                            echo "chemistry,ARC-v1" >> "${library_output}"
                        fi
                        echo "" >> "${library_output}"
                        echo "[vdj]" >> "${library_output}"
                        echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-IMGT-VDJ-2024" >> "${library_output}"
                        # Add options if there are VDJ-specific options specified by --vdj-options
                        if [ "${vdj_options}" != "NA" ]; then
                            IFS=',' read -ra values <<< "${vdj_options}"
                            for value in "${values[@]}"; do
                                echo "${value}" >> "${library_output}"
                            done
                        fi
                    # Is it a mouse run?
                    elif [[ "${species}" =~ ^(Mouse|mouse|Mm|mm)$ ]]; then
                        echo "[gene-expression]" >> "${library_output}"
                        echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-hardmasked-optimised-arc" >> "${library_output}"
#                       echo "probe-set,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/frp-probes/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" >> "${library_output}"
                        echo "create-bam,true" >> "${library_output}"
                        # Add options if there are gene expression-specific options specified by --gene-expression-options
                            if [ -n "${gene_expression_options}" ] && [ "${gene_expression_options}" != "NA" ]; then
                                IFS=';' read -ra values <<< "${gene_expression_options}"
                                for value in "${values[@]}"; do
                                    echo "${value}" >> "${library_output}"
                                done
                            fi
                        # If this is a DOGMA-seq or MULTIOME run, to specify chemistry
                        if [ "${assay}" == "DOGMA" ] || [ "${assay}" == "MULTIOME" ]; then
                            echo "chemistry,ARC-v1" >> "${library_output}"
                        fi
                        echo "" >> "${library_output}"
                        echo "[vdj]" >> "${library_output}"
                        echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-IMGT-VDJ-2024" >> "${library_output}"
                        # Add options if there are VDJ-specific options specified by --vdj-options
                        if [ -n "${vdj_options}" ] && [ "${vdj_options}" != "NA" ]; then
                            IFS=',' read -ra values <<< "${vdj_options}"
                            for value in "${values[@]}"; do
                                echo "${value}" >> "${library_output}"
                            done
                        fi
                    fi
                    # Is there paired ADT data for this sample?
                    if [ "${adt_file}" != "NA" ]; then
                        echo "" >> "${library_output}"
                        echo "[feature]" >> "${library_output}"
                        echo "reference,$project_scripts/ADT_files/${adt_file}.csv" >> "${library_output}"
                        # Add options if there are ADT/HTO-specific options specified by --adt-options
                        if [ "${adt_options}" != "" ]; then
                            IFS=',' read -ra values <<< "${adt_options}"
                            for value in "${values[@]}"; do
                                echo "${value}" >> "${library_output}"
                            done
                        fi
                    fi
                    echo "Writing ${modality} for ${library}"
                    echo "" >> "${library_output}"
                    echo "[libraries]" >> "${library_output}"
                    echo "fastq_id,fastqs,feature_types" >> "${library_output}"
                fi
            elif [[ "${modality}" == 'HTO' || "${modality}" == 'ADT' || "${modality}" == 'VDJ-T' || "${modality}" == 'VDJ-B' || "${modality}" == 'CRISPR' ]]; then
                if [[ -f ${library_output} ]]; then
                    :
                else
                    echo -e "\033[0;31mERROR:\033[0m Please ensure that in the metadata file, GEX libraries for all samples are first, before ADT/HTO/VDJ-T/CRISPR"
                    exit 1
                fi
            fi
            # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
            declare -A unique_lines
            # Recursively search for FASTQ files in the project_id FASTQ folder
            for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
                matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
                echo $matching_fastq_files
                for fastq_file in "${matching_fastq_files[@]}"; do
                    # Extract the directory containing the FASTQ file
                    directory=$(dirname "${fastq_file}")
                    # Extract the modified name from the FASTQ file
                    fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                    # Create a unique identifier the FASTQ file
                    line_identifier="${fastq_name},${directory},${full_modality}"
                    # Check if the line has already been added
                    if [ ! -v unique_lines["${line_identifier}"] ]; then
                        unique_lines["${line_identifier}"]=1
                        echo "${fastq_name},${directory},${full_modality}" >> "${library_output}"
                        echo "Writing ${fastq_name},${directory},${full_modality} to ${library}"
                    fi
                done
            done

        # If this run is of ATAC samples, alongside ADT/HTO
        elif [[ ${run_type} == 'ATAC' && ${mode} == "GEX" ]]; then
            if [[ (${modality} == 'ADT' || ${modality} == 'HTO') && ${assay} != 'ASAP' ]]; then
                if [[ -f ${library_output} ]]; then
                    :
                else
                    echo -e "\033[0;31mERROR:\033[0m If you're trying to combine DOGMA ADT/HTO to a DOGMA GEX, please make sure the output directory is of the run containing the GEX FASTQ files"
                    echo -e "\033[0;31mERROR:\033[0m The reference part of the csv file needs to be initialised"
                    exit 1
                fi
                # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
                declare -A unique_lines
                # Recursively search for FASTQ files in the project_id FASTQ folder
                for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
                    matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
                    for fastq_file in "${matching_fastq_files[@]}"; do
                        # Extract the directory containing the FASTQ file
                        directory=$(dirname "${fastq_file}")
                        # Extract the modified name from the FASTQ file
                        fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                        # Create a unique identifier the FASTQ file
                        line_identifier="${fastq_name},${directory},${full_modality}"
                        # Check if the line has already been added
                        if [ ! -v unique_lines["${line_identifier}"] ]; then
                            unique_lines["${line_identifier}"]=1
                            echo "${fastq_name},${directory},${full_modality}" >> "${library_output}"
                            echo "Writing ${fastq_name},${directory},${full_modality} to ${library}"
                        fi
                    done
                done
            fi
        elif [[ ${run_type} == 'ATAC' && ${mode} == "ATAC" ]]; then
            if [[ (${modality} == 'ADT' || ${modality} == 'HTO') && ${assay} == 'DOGMA' ]]; then
                :
            elif [[ (${modality} == 'ADT' || ${modality} == 'HTO') && ${assay} == 'ASAP' ]]; then
                library_output=${output_project_libraries}/${library}_ADT.csv
                # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
                declare -A unique_lines
                # Recursively search for FASTQ files in the project_id FASTQ folder
                for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
                    matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
                    for fastq_file in "${matching_fastq_files[@]}"; do
                        # Extract the directory containing the FASTQ file
                        directory=$(dirname "${fastq_file}")
                        # Extract the modified name from the FASTQ file
                        fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                        # Create a unique identifier the FASTQ file
                        line_identifier="${fastq_name},${directory},${full_modality}"
                        # Check if the line has already been added
                        if [ ! -v unique_lines["${line_identifier}"] ]; then
                            unique_lines["${line_identifier}"]=1
                            echo "${fastq_name},${directory}" >> "${library_output}"
                            echo "Writing ${fastq_name},${directory} to ${library}"
                        fi
                    done
                done
            elif [[ ${modality} == 'ATAC' ]]; then
                # Initialize an associative array, as the script works by checking for wildcard sample name of FASTQ files and only one sample will be added per FASTQ group
                declare -A unique_lines
                # Recursively search for FASTQ files in the project_id FASTQ folder
                for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
                    matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
                    for fastq_file in "${matching_fastq_files[@]}"; do
                        # Extract the directory containing the FASTQ file
                        directory=$(dirname "${fastq_file}")
                        # Extract the modified name from the FASTQ file
                        fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
                        # Create a unique identifier the FASTQ file
                        line_identifier="${fastq_name},${directory},${full_modality}"
                        # Check if the line has already been added
                        if [ ! -v unique_lines["${line_identifier}"] ]; then
                            unique_lines["${line_identifier}"]=1
                            echo "${fastq_name},${directory}" >> "${library_output}"
                            echo "Writing ${fastq_name},${directory} to ${library}.csv"
                        fi
                    done
                done
            else
                echo -e "\033[0;31mERROR:\033[0m Cannot determine modality for this ATAC run. Are you sure the only modalities are either ATAC, ADT, or HTO?"
                echo -e "\033[0;31mERROR:\033[0m Library: ${library}, modality: ${modality}"
                exit 1
            fi
        fi
    done < "${metadata_file}"
done

# Ask the user if they want to submit the indices for FASTQ generation
echo "Would you like to proceed to counting? (Y/N)"
read -r choice
while [[ ! ${choice} =~ ^[YyNn]$ ]]; do
    echo "Invalid input. Please enter y or n"
    read -r choice
done

# Process choices
if [ "${choice}" = "Y" ] || [ "${choice}" = "y" ]; then
    echo "Submitting: bash ${oscar_dir}/04_count.sh --project-id ${output_project_id}"
    bash ${oscar_dir}/04_count.sh --project-id ${output_project_id}
elif [ "${choice}" = "N" ] || [ "${choice}" = "n" ]; then
    :
else
    echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
fi
