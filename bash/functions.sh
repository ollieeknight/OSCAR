#!/bin/bash

# Function for formatted logging
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "[$timestamp] $1"
}

# Function to check last command status and exit if failed
check_status() {
    if [ $? -eq 0 ]; then
        log "✓ SUCCESS: $1"
    else
        log "❌ ERROR: $1 failed"
        exit 1
    fi
}

# Function to check if project_id is defined
check_project_id() {
    if [ -z "${project_id}" ]; then
        echo -e "\033[0;31mERROR:\033[0m project-id is not defined. Please provide --project-id"
        exit 1
    fi
}

check_folder_exists() {
    folder_name="$1"
    if [ ! -d "$folder_name" ]; then
        echo -e "\033[0;31mERROR:\033[0m Folder '$folder_name' does not exist. Please provide a valid folder name."
        exit 1
    fi
}

# Function to check run_type based on read length from RunInfo.xml
check_run_type() {
    project_id="$1"
    dir_prefix="$2"

    read_length=$(awk -F '"' '/<Read Number="1"/ {print $4}' "${dir_prefix}/${project_id}/${project_id}_bcl/RunInfo.xml")

    if [ -z "${read_length}" ]; then
        echo -e "\033[0;31mERROR:\033[0m Unable to find read length. Please check ${dir_prefix}/${project_id}_bcl/RunInfo.xml"
        exit 1
    fi

    if [ "${read_length}" -gt 45 ]; then
        run_type="ATAC"
    elif [ "${read_length}" -lt 45 ]; then
        run_type="GEX"
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine run type, please check ${dir_prefix}/${project_id}_bcl/RunInfo.xml"
        exit 1
    fi

    echo "${run_type}"
}

# Function to check if metadata file exists
check_metadata_file() {
    metadata_file="$1"
    if [[ ! -f "${metadata_file}" ]]; then
        echo -e "\033[0;31mERROR:\033[0m Metadata file not found under ${metadata_file}"
        exit 1
    fi
}

check_and_pull_oscar_containers() {

    container_count="${TMPDIR}/OSCAR/oscar-count_latest.sif"

    if [ ! -f "${container_count}" ]; then
        echo "oscar-count_latest.sif singularity file not found, pulling..."
        mkdir -p "${TMPDIR}/OSCAR"
        apptainer pull --dir "${TMPDIR}/OSCAR" library://romagnanilab/oscar/oscar-count:latest
    fi

    container_qc="${TMPDIR}/OSCAR/oscar-qc_latest.sif"

    if [ ! -f "${container_qc}" ]; then
        echo "oscar-qc_latest.sif singularity file not found, pulling..."
        mkdir -p "${TMPDIR}/OSCAR"
        apptainer pull --dir "${TMPDIR}/OSCAR" library://romagnanilab/oscar/oscar-qc:latest
        echo "All images are present under ${TMPDIR}/OSCAR/"
        touch ${TMPDIR}/OSCAR/oscar-count_latest.sif
        touch ${TMPDIR}/OSCAR/oscar-qc_latest.sif
    else
        echo "All images are present under ${TMPDIR}/OSCAR/"
        touch ${TMPDIR}/OSCAR/oscar-count_latest.sif
        touch ${TMPDIR}/OSCAR/oscar-qc_latest.sif
    fi
}

check_base_masks_step1() {
    local xml_file="${project_dir}/${project_id}_bcl/RunInfo.xml"

    if [ ! -f "${xml_file}" ]; then
        echo "Sequencing run RunInfo.xml file not found under ${project_dir}/${project_id}_bcl/. This is required to determine base masks."
        exit 1
    fi

    if grep -q '<Reads>' "${xml_file}" && grep -q '</Reads>' "${xml_file}"; then
        num_reads=$(grep -o '<Read Number="' "${xml_file}" | wc -l)
        if [ "${num_reads}" -eq 3 ]; then
            reads=3
        elif [ "${num_reads}" -eq 4 ]; then
            reads=4
        else
            echo -e "\033[0;31mERROR:\033[0m RunInfo.xml contains unexpected reads, please check the file"
            exit 1
        fi
    else
        echo -e "\033[0;31mERROR:\033[0m RunInfo.xml is missing expected <Reads> tags"
        exit 1
    fi
}

check_base_masks_step2() {
    if [[ ${reads} == 3 ]]; then
        base_mask_SI_3prime_GEX='Y28n*,I8n*,Y90n*'
        base_mask_DI_3prime_GEX='Y28n*,I8n*,Y90n*'
        base_mask_SI_3prime_ADT='Y28n*,I8n*,Y90n*'
        base_mask_DOGMA_ADT='Y28n*,I8n*,Y90n*'
    elif [[ ${reads} == 4 ]]; then
        base_mask_SI_3prime_GEX='Y28n*,I8n*,N*,Y90n*'
        base_mask_SI_5prime_GEX='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DI_3prime_GEX='Y28n*,I8n*,I8n*,Y90n*'
        base_mask_DI_5prime_GEX='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_SI_3prime_ADT='Y28n*,I8n*,N*,Y90n*'
        base_mask_DI_5prime_ADT='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DOGMA_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DOGMA_ATAC='Y100n*,I8n*,Y24n*,Y100n*'
        base_mask_DOGMA_ADT='Y28n*,I8n*,N*,Y90n*'
        base_mask_ATAC_ATAC='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_ASAP_ATAC='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_ASAP_ADT='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_ASAP_GENO='Y100n*,I8n*,Y16n*,Y100n*'
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine number of reads, check RunInfo.xml file"
        exit 1
    fi
}

check_base_masks_step3() {
    local file="$1"
    local run_type="$2"

    # Default values
    local cellranger_command=""
    local index_type=""
    local filter_option=""
    local base_mask=""

    # Logic for determining the parameters
    if [[ ${file} == CITE* ]]; then
        cellranger_command='cellranger mkfastq'
        if [[ ${file} == *_SI_* ]]; then
            index_type='SI'
            filter_option='--filter-single-index'
            if [[ ${file} == *_GEX* ]]; then
                base_mask=$base_mask_SI_3prime_GEX
            elif [[ ${file} == *_ADT* ]] || [[ ${file} == *_HTO* ]]; then
                base_mask=$base_mask_SI_3prime_ADT
            fi
        elif [[ ${file} == *_DI_* ]]; then
            index_type='DI'
            filter_option='--filter-dual-index'
            if [[ ${file} == *_3prime* ]]; then
                if [[ ${file} == *_GEX* ]]; then
                    base_mask=$base_mask_DI_3prime_GEX
                elif [[ ${file} == *_ADT* ]] || [[ ${file} == *_HTO* ]]; then
                    base_mask=$base_mask_DI_3prime_ADT
                fi
            elif [[ ${file} == *_5prime* ]]; then
                if [[ ${file} == *_GEX* ]] || [[ ${file} == *_VDJ-* ]]; then
                    base_mask=$base_mask_DI_5prime_GEX
                elif [[ ${file} == *_ADT* ]] || [[ ${file} == *_HTO* ]]; then
                    base_mask=$base_mask_DI_5prime_ADT
                fi
            fi
        fi
    elif [[ ${file} == GEX* ]]; then
        cellranger_command='cellranger mkfastq'
        if [[ ${file} == *_SI* ]]; then
            index_type='SI'
            filter_option='--filter-single-index'
            base_mask=$base_mask_SI_3prime_GEX
        elif [[ ${file} == *_DI* ]]; then
            index_type='DI'
            filter_option='--filter-dual-index'
            if [[ ${file} == *_3prime* ]]; then
                base_mask=$base_mask_DI_3prime_GEX
            elif [[ ${file} == *_5prime* ]]; then
                base_mask=$base_mask_DI_5prime_GEX
            fi
        fi
    elif [[ ${file} == DOGMA* ]]; then
        index_type='DI'
        filter_option='--filter-dual-index'
        if [[ ${file} == *_GEX* ]]; then
            cellranger_command='cellranger mkfastq'
            base_mask=$base_mask_DOGMA_GEX
        elif [[ ${file} == *_ADT* ]] || [[ ${file} == *_HTO* ]]; then
            base_mask=$base_mask_DOGMA_ADT
            cellranger_command='cellranger mkfastq'
            [[ ${run_type} == 'ATAC' ]] && cellranger_command='cellranger-atac mkfastq'
        elif [[ ${file} == *_ATAC* ]]; then
            cellranger_command='cellranger-atac mkfastq'
            base_mask=$base_mask_DOGMA_ATAC
        fi
    elif [[ ${file} == ASAP* ]]; then
        cellranger_command='cellranger-atac mkfastq'
        index_type='DI'
        filter_option='--filter-dual-index'
        if [[ ${file} == *_ATAC* ]]; then
            base_mask=$base_mask_ASAP_ATAC
        elif [[ ${file} == *_ADT* ]] || [[ ${file} == *_HTO* ]]; then
            base_mask=$base_mask_ASAP_ADT
        elif [[ ${file} == *_GENO* ]]; then
            base_mask=$base_mask_ASAP_GENO
        fi
    elif [[ ${file} == ATAC* ]]; then
        cellranger_command='cellranger-atac mkfastq'
        index_type='DI'
        filter_option='--filter-dual-index'
        base_mask=$base_mask_ASAP_ATAC
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine base mask for ${file}, please check path"
        exit 1
    fi

    # Export the variables if needed
  echo "${cellranger_command// /.}" "${index_type// /.}" "${filter_option// /.}" "${base_mask// /.}"
}

validate_mode() {
    local mode=$1

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
}

print_options() {
    local options=("$@")
    local option_name=$1
    shift
    local option_values=("$@")

    if [ ${#option_values[@]} -gt 0 ]; then
        echo "${option_name} options set as:"
        for option in "${option_values[@]}"; do
            echo "${option}"
        done
    else
        echo "No non-default options set for ${option_name}"
    fi
}

determine_full_modality() {
    local modality=$1
    local library=$2
    local full_modality=""

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
        return 1
    fi

    echo "${full_modality}"
}

write_human_reference() {
    echo -e "\033[0;33mWriting human reference files for ${library}\033[0m"
    echo "[gene-expression]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-hardmasked-optimised-arc" >> "${library_output}"
    echo "create-bam,true" >> "${library_output}"
    if [ -n "${gene_expression_options}" ] && [ "${gene_expression_options}" != "NA" ]; then
        IFS=';' read -ra values <<< "${gene_expression_options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
    if [ "${assay}" == "DOGMA" ] || [ "${assay}" == "MULTIOME" ]; then
        echo "chemistry,ARC-v1" >> "${library_output}"
    fi
    echo "" >> "${library_output}"
    echo "[vdj]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-IMGT-VDJ-2024" >> "${library_output}"
    if [ "${vdj_options}" != "NA" ]; then
        IFS=',' read -ra values <<< "${vdj_options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
}

write_mouse_reference() {
    echo -e "\033[0;33mWriting mouse reference files for ${library}\033[0m"
    echo "[gene-expression]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-hardmasked-optimised-arc" >> "${library_output}"
    echo "create-bam,true" >> "${library_output}"
    if [ -n "${gene_expression_options}" ] && [ "${gene_expression_options}" != "NA" ]; then
        IFS=';' read -ra values <<< "${gene_expression_options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
    if [ "${assay}" == "DOGMA" ] || [ "${assay}" == "MULTIOME" ]; then
        echo "chemistry,ARC-v1" >> "${library_output}"
    fi
    echo "" >> "${library_output}"
    echo "[vdj]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-IMGT-VDJ-2024" >> "${library_output}"
    if [ -n "${vdj_options}" ] && [ "${vdj_options}" != "NA" ]; then
        IFS=',' read -ra values <<< "${vdj_options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
}

write_adt_data() {
    echo "" >> "${library_output}"
    echo "[feature]" >> "${library_output}"
    echo "reference,$project_scripts/ADT_files/${adt_file}.csv" >> "${library_output}"
    if [ "${adt_options}" != "" ]; then
        IFS=',' read -ra values <<< "${adt_options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
}

write_fastq_files() {
    declare -A unique_lines
    for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
        matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
        for fastq_file in "${matching_fastq_files[@]}"; do
            directory=$(dirname "${fastq_file}")
            fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//' | sed -E 's/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
            
            # Define line_identifier and output suffix based on conditions
            if [[ "${modality}" == "ADT" && "${assay}" == "ASAP" ]]; then
                line_identifier="${fastq_name},${directory}"
                suffix="_ADT"
            elif [[ "${modality}" == "ATAC" || "${assay}" == "ASAP" ]]; then
                line_identifier="${fastq_name},${directory},${full_modality}"
                suffix="_ATAC"
            else
                line_identifier="${fastq_name},${directory},${full_modality}"
                suffix=""
            fi

            # Construct the final output file name
            output_file="${library_output%.csv}"
            if [[ ! "${output_file}" =~ ${suffix}$ ]]; then
                output_file="${output_file}${suffix}"
            fi
            if [[ ! "${output_file}" =~ \.csv$ ]]; then
                output_file="${output_file}.csv"
            fi

            if [ ! -v unique_lines["${line_identifier}"] ]; then
                unique_lines["${line_identifier}"]=1
                echo "${line_identifier}" >> "${output_file}"
                echo -e "\033[0;33mWriting ${line_identifier} to\033[0m"
                echo "${output_file}"
            fi
        done
    done
}

handle_gex_mode() {
    if [[ "${modality}" == 'GEX' ]]; then
        if [ ! -f "${library_output}" ]; then
            if [[ "${species}" =~ ^(Human|human|Hs|hs)$ ]]; then
                write_human_reference
            elif [[ "${species}" =~ ^(Mouse|mouse|Mm|mm)$ ]]; then
                write_mouse_reference
            fi
            if [ "${adt_file}" != "NA" ]; then
                write_adt_data
            fi
            echo -e "\033[0;33mWriting ${modality} for ${library}\033[0m"
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
    write_fastq_files
}

handle_atac_mode() {
    if [[ (${modality} == 'ADT' || ${modality} == 'HTO') && ${assay} != 'ASAP' ]]; then
        if [[ -f ${library_output} ]]; then
            :
        else
            echo -e "\033[0;31mERROR:\033[0m If you're trying to combine DOGMA ADT/HTO to a DOGMA GEX, please make sure the output directory is of the run containing the GEX FASTQ files"
            echo -e "\033[0;31mERROR:\033[0m The reference part of the csv file needs to be initialised"
            exit 1
        fi
        write_fastq_files
    elif [[ (${modality} == 'ADT' || ${modality} == 'HTO') && ${assay} == 'ASAP' ]]; then
        library_output=${output_project_libraries}/${library}_ADT.csv
        write_fastq_files
    elif [[ ${modality} == 'ATAC' ]]; then
        write_fastq_files
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine modality for this ATAC run. Are you sure the only modalities are either ATAC, ADT, or HTO?"
        echo -e "\033[0;31mERROR:\033[0m Library: ${library}, modality: ${modality}"
        exit 1
    fi
}

count_read_csv() {
    local library_folder=$1
    local library=$2
    local fastq_names=""
    local fastq_dirs=""

    while IFS= read -r line; do
        if [[ $line == *ATAC* ]]; then
            IFS=',' read -r fastq_name fastq_dir <<< "$line"
            fastq_names="${fastq_names:+$fastq_names,}${fastq_name}"
            fastq_dirs="${fastq_dirs:+$fastq_dirs,}${fastq_dir}"
        fi
    done < "${library_folder}/${library}.csv"

    echo "$fastq_names" "$fastq_dirs"
}

count_check_dogma() {
    local library_folder=$1
    local library=$2
    local extra_arguments=""

    if grep -q '.*DOGMA.*' "${library_folder}/${library}.csv"; then
        cat "${library_folder}/${library}.csv"
        extra_arguments="--chemistry ARC-v1"
        echo "Adding $extra_arguments as it is a DOGMA/MULTIOME run"
    fi

    echo "$extra_arguments"
}

count_read_metadata() {
    local metadata_file=$1
    local library=$2
    local ADT_file=""

    while IFS=',' read -r -a fields; do
        assay="${fields[0]}"
        experiment_id="${fields[1]}"
        historical_number="${fields[2]}"
        replicate="${fields[3]}"
        expected_library="${assay}_${experiment_id}_exp${historical_number}_lib${replicate}_ATAC"

        if [ "$expected_library" == "$library" ]; then
            ADT_file="${fields[10]}"  # Corrected index for ADT file
            break
        fi
    done < "$metadata_file"

    echo "$ADT_file"
}

count_read_adt_csv() {
    local library_folder=$1
    local library=$2
    local library_csv="${library_folder}/${library}"
    local fastq_dirs=''
    local fastq_libraries=''

    if [[ ! -f "${library_csv}" ]]; then
        echo "ERROR: ${library_csv} not found."
        exit 1
    fi

    while IFS= read -r line; do
        IFS=',' read -r -a parts <<< "$line"
        fastq_library="${parts[0]}"
        fastq_dir="${parts[1]}"
        fastq_libraries="${fastq_libraries},${fastq_library}"
        fastq_dirs="${fastq_dirs},${fastq_dir}"
    done < "${library_csv}"

    fastq_dirs="${fastq_dirs:1}"
    fastq_libraries="${fastq_libraries:1}"

    echo "$fastq_dirs" "$fastq_libraries"
}

extract_variables() {
    local library="$1"
    local assay remainder experiment_id historical_number replicate

    assay="${library%%_*}"
    remainder="${library#*_}"
    experiment_id="${remainder%%_exp*}"
    remainder="${remainder#*_}"
    historical_number="${remainder%%_lib*}"
    replicate="${remainder#*_lib}"

    echo "$assay" "$experiment_id" "$historical_number" "$replicate"
}

search_metadata() {
    local library="$1"
    local assay="$2"
    local experiment_id="$3"
    local historical_number="$4"
    local replicate="$5"
    local project_ids=("${!6}")
    local prefix="$7"
    local n_donors ADT_file

    # Loop through each project_id
    for project_id in "${project_ids[@]}"; do
        metadata_file="${prefix}/${project_id}/${project_id}_scripts/metadata/metadata.csv"

        # Check if the metadata file exists
        if [ -f "$metadata_file" ]; then
            echo "Searching $project_id for ${library}"
            # Read metadata from the CSV file line by line
            while IFS=, read -r -a fields; do
                # Check if all individual fields match the criteria
                if [[ "${fields[0]}" == "$assay" && "${fields[1]}" == "$experiment_id" && "${fields[2]}" == "$historical_number" && "${fields[3]}" == "$replicate" ]]; then
                    n_donors="${fields[9]}"
                    ADT_file="${fields[10]}"
                    break  # Stop searching once a match is found
                fi
            done < "$metadata_file"
        else
            echo -e "\033[0;31mERROR:\033[0m Metadata file not found for project_id: $project_id"
            exit 1
        fi
    done

    echo "$n_donors" "$ADT_file"
}