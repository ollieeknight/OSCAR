#!/bin/bash

# Function to check if project_id is defined
check_project_id() {
    if [ -z "${project-id}" ]; then
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
    container="${TMPDIR}/OSCAR/oscar-count_latest.sif"
    
    if [ ! -f "${container}" ]; then
        echo "oscar-count_latest.sif singularity file not found, pulling..."
        mkdir -p "${TMPDIR}/OSCAR"
        apptainer pull --dir "${TMPDIR}/OSCAR" library://romagnanilab/default/oscar-count:latest
    fi

        # container="${TMPDIR}/OSCAR/oscar-qc_latest.sif"
    
    # if [ ! -f "${container}" ]; then
    #     echo "oscar-qc_latest.sif singularity file not found, pulling..."
    #     mkdir -p "${TMPDIR}/OSCAR"
    #     apptainer pull --dir "${TMPDIR}/OSCAR" library://romagnanilab/default/oscar-qc:latest
    # fi
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

# Function to check and print options for a given array
print_options() {
    local array_name="$1"
    local array_ref=("${!2}")
    local label="$3"

    if [ ${#array_ref[@]} -gt 0 ]; then
        echo "${label} set as:"
        for option in "${array_ref[@]}"; do
            echo "${option}"
        done
    else
        echo "No options set for ${label,,}"  # Convert label to lowercase
    fi
}

process_library_files() {
    local library_type="$1"
    local library_output="$2"

    if [[ "${run_type}" == "${library_type}" && "${mode}" == "${library_type}" ]]; then
        if [[ "${modality}" == 'GEX' ]]; then
            handle_gex "${library_output}"
        elif [[ "${modality}" =~ ^(ADT|HTO|VDJ-T|VDJ-B|CRISPR)$ ]]; then
            handle_other "${library_output}"
        fi
    fi
}

handle_gex() {
    local library_output="$1"
    # Check for existing output file
    if [[ -f "${library_output}" ]]; then
        return
    fi

    echo "Writing ${modality} for ${library}"

    if [[ "${species}" =~ ^(Human|human|Hs|hs)$ ]]; then
        write_human_gex "${library_output}"
    elif [[ "${species}" =~ ^(Mouse|mouse|Mm|mm)$ ]]; then
        write_mouse_gex "${library_output}"
    fi

    # Handle ADT data if specified
    if [[ "${adt_file}" != "NA" ]]; then
        echo "" >> "${library_output}"
        echo "[feature]" >> "${library_output}"
        echo "reference,$project_scripts/ADT_files/${adt_file}.csv" >> "${library_output}"
        if [[ -n "${adt_options}" ]]; then
            IFS=',' read -ra values <<< "${adt_options}"
            for value in "${values[@]}"; do
                echo "${value}" >> "${library_output}"
            done
        fi
    fi

    echo "" >> "${library_output}"
    echo "[libraries]" >> "${library_output}"
    echo "fastq_id,fastqs,feature_types" >> "${library_output}"

    # Process fastq files
    write_fastq_files "${library_output}"
}

write_human_gex() {
    local library_output="$1"
    echo "[gene-expression]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-hardmasked-optimised-arc" >> "${library_output}"
    echo "create-bam,true" >> "${library_output}"
    add_options "${library_output}" "${gene_expression_options}"

    if [[ "${assay}" == "DOGMA" || "${assay}" == "MULTIOME" ]]; then
        echo "chemistry,ARC-v1" >> "${library_output}"
    fi

    echo "" >> "${library_output}"
    echo "[vdj]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-IMGT-VDJ-2024" >> "${library_output}"
    add_options "${library_output}" "${vdj_options}"
}

write_mouse_gex() {
    local library_output="$1"
    echo "[gene-expression]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-hardmasked-optimised-arc" >> "${library_output}"
    echo "create-bam,true" >> "${library_output}"
    add_options "${library_output}" "${gene_expression_options}"

    if [[ "${assay}" == "DOGMA" || "${assay}" == "MULTIOME" ]]; then
        echo "chemistry,ARC-v1" >> "${library_output}"
    fi

    echo "" >> "${library_output}"
    echo "[vdj]" >> "${library_output}"
    echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/mm/GRCm38-IMGT-VDJ-2024" >> "${library_output}"
    add_options "${library_output}" "${vdj_options}"
}

add_options() {
    local library_output="$1"
    local options="$2"
    if [[ -n "${options}" && "${options}" != "NA" ]]; then
        IFS=';' read -ra values <<< "${options}"
        for value in "${values[@]}"; do
            echo "${value}" >> "${library_output}"
        done
    fi
}

write_fastq_files() {
    local library_output="$1"
    declare -A unique_lines
    for folder in "${project_dir}/${project_id}_fastq"/*/outs; do
        matching_fastq_files=($(find "${folder}" -type f -name "${library}*${modality}*" | sort -u))
        for fastq_file in "${matching_fastq_files[@]}"; do
            directory=$(dirname "${fastq_file}")
            fastq_name=$(basename "${fastq_file}" | sed -E 's/\.fastq\.gz$//; s/(_S[0-9]+)?(_[SL][0-9]+_[IR][0-9]+_[0-9]+)*$//')
            line_identifier="${fastq_name},${directory},${full_modality}"

            if [ ! -v unique_lines["${line_identifier}"] ]; then
                unique_lines["${line_identifier}"]=1
                echo "${fastq_name},${directory},${full_modality}" >> "${library_output}"
                echo "Writing ${fastq_name},${directory},${full_modality} to ${library_output}"
            fi
        done
    done
}
