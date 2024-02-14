#!/bin/bash

# usage is: bash 05_quality_control.sh -project-id project_ids

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
# │   ├── FASTQ_2
# │   ├── ...
# │   ├── FASTQ_n
# └── {project_id}_scripts/
#     ├── ADT_files/
#     ├── indices/
#     ├── libraries/
#     └── metadata/ # METADATA MUST BE IN THIS FOLDER!

# Define default values
OSCAR_script_dir=$(dirname "${BASH_SOURCE[0]}")
prefix="$HOME/scratch/ngs"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --project-id)
        project_ids_in="$2"
        shift 2
        ;;
        --prefix)
        prefix="$2"
        shift 2
        ;;
        *)
    esac
done

IFS=',' read -r -a project_ids <<< "$project_ids_in"

output_project_id=${project_ids[0]}

# Check if project_id is empty
if [ -z "$output_project_id" ]; then
    echo "Please provide a project_id using the --project-id option."
    exit 1
fi

# Define the project directory
output_project_dir=${prefix}/$output_project_id
outs="${prefix}/${output_project_id}/${output_project_id}_outs"

# Check if the output_project_id folder exists
if [ ! -d "$outs" ]; then
    echo "Error: Outs folder ($outs) does not exist. Make sure cellranger has run"
    exit 1
fi
container="$TMPDIR/oscar-qc_latest.sif"

# Check if the container exists in TMPDIR
if [ ! -f "$container" ]; then
    # Check if the container exists in the alternative directory
    if [ -f "$HOME/group/work/bin/OSCAR/singularity/oscar-qc_latest.sif" ]; then
        # Copy the container to TMPDIR
        cp "$HOME/group/work/bin/OSCAR/singularity/oscar-qc_latest.sif" "$TMPDIR/"
    else
        echo "oscar-qc_latest.sif singularity file not found, pulling..."
        apptainer pull --dir $TMPDIR library://romagnanilab/default/oscar-qc:latest
    fi
fi

# Check if the container is now available in TMPDIR
if [ -f "$container" ]; then
    echo "oscar-qc_latest.sif Singularity file found in $TMPDIR."
else
    echo "Failed to find oscar-qc_latest.sif in $TMPDIR."
fi

libraries=($(find ${prefix}/$output_project_id/${output_project_id}_outs/ -maxdepth 1 -mindepth 1 -type d -not -name 'logs' -exec basename {} \;))

for library in "${libraries[@]}"; do
    echo ""
    echo "----"
    echo ""
    # Extract variables from the library name
    assay="${library%%_*}"
    remainder="${library#*_}"
    experiment_id="${remainder%%_exp*}"
    remainder="${remainder#*_}"
    remainder="${remainder#*exp}"
    historical_number="${remainder%%_lib*}"
    remainder="${remainder#*lib}"
    replicate="${remainder}"

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
            echo "Error: Metadata file not found for project_id: $project_id"
            exit 1
        fi
    done

    feature_matrix_path=$(find "${outs}/${library}/" -type f -name "raw_feature_bc_matrix.h5" -print -quit)
    peak_matrix_path=$(find "${outs}/${library}/" -type f -name "raw_peak_bc_matrix.h5" -print -quit)

    if [ -n "$feature_matrix_path" ]; then
        read -p "Would you like to submit ambient RNA removal with cellbender for ${library}? (Y/N)" perform_function
        # Convert input to uppercase for case-insensitive comparison
        perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
        if [ "$perform_function" != "Y" ]; then
            echo "Skipping cellbender for ${library}"
        else
            mkdir -p ${outs}/${library}/cellbender/
job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_cellbender
#SBATCH --output --output $outs/logs/${library}_cellbender.out
#SBATCH --error $outs/logs/${library}_cellbender.out
#SBATCH --ntasks 1
#SBATCH --partition "gpu"
#SBATCH --gres gpu:1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64000
#SBATCH --time 5:00:00
cd ${outs}/${library}
apptainer run --nv -B /fast ${container} cellbender remove-background --cuda --input ${feature_matrix_path} --output ${outs}/${library}/cellbender/output.h5
rm ckpt.tar.gz
EOF
        )
        job_id=$(echo "$job_id" | awk '{print $4}')
        echo ""
        fi
        if [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            echo "Skipping genotyping for ${library}, as this is either a mouse run, or only contains 1 donor"
            job_id=""
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" != "" ]]; then
             read -p "Would you like to genotype ${library}? (Y/N)" perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping for ${library}"
            else
                echo "Submitting vireo genotyping for ${library}"
                mkdir -p ${outs}/${library}/vireo/
sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_vireo
#SBATCH --output $outs/logs/${library}_vireo.out
#SBATCH --error $outs/logs/${library}_vireo.out
#SBATCH --ntasks=32
#SBATCH --mem=32000
#SBATCH --time=6:00:00
# The following line ensures that this job runs after the previous job with ID $job_id
#SBATCH --dependency=afterok:$job_id
num_cores=\$(nproc)
cd ${outs}/${library}
apptainer run -B /fast ${container} cellsnp-lite -s ${outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam -b ${outs}/${library}/cellbender/output_cell_barcodes.csv -O ${outs}/${library}/vireo -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores
apptainer run -B /fast ${container} vireo -c ${outs}/${library}/vireo -o ${outs}/${library}/vireo -N $n_donors -p \$num_cores
EOF
                job_id=""
            fi
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" == "" ]]; then
             read -p "Would you like to genotype ${library}? (Y/N)" perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping for ${library}"
            else
                echo "Submitting vireo genotyping for ${library}"
                mkdir -p ${outs}/${library}/vireo/
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_vireo
#SBATCH --output $outs/logs/${library}_vireo.out
#SBATCH --error $outs/logs/${library}_vireo.out
#SBATCH --ntasks=32
#SBATCH --mem=32000
#SBATCH --time=6:00:00
num_cores=\$(nproc)
cd ${outs}/${library}
apptainer run -B /fast ${container} cellsnp-lite -s ${outs}/${library}/outs/per_sample_outs/${library}/count/sample_alignments.bam -b ${outs}/${library}/cellbender/output_cell_barcodes.csv -O ${outs}/${library}/vireo -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores
apptainer run -B /fast ${container} vireo -c ${outs}/${library}/vireo -o ${outs}/${library}/vireo -N $n_donors -p \$num_cores
EOF
                job_id=""
            fi
        else
            echo "Cannot determine the number of donors for ${library}"
            exit 1
        fi
    elif [ -n "$peak_matrix_path" ]; then
        if [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' ]]; then
            read -p "Would you like to genotype ${library}? (Y/N): " perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')

            # Check if the input is 'N' or 'n'
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping"
            else
                mkdir -p ${outs}/${library}/AMULET
                mkdir -p ${outs}/${library}/vireo/
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_QC
#SBATCH --output $outs/logs/${library}_genotyping.out
#SBATCH --error $outs/logs/${library}_genotyping.out
#SBATCH --ntasks=32
#SBATCH --mem=68000
#SBATCH --time=8:00:00
num_cores=\$(nproc)
cd ${outs}/${library}
echo "Starting mgatk mtDNA genotyping"
echo ""
apptainer run -B /fast,/usr ${container} mgatk tenx -i ${outs}/${library}/outs/possorted_bam.bam -n output -o ${outs}/${library}/mgatk -c 8 -bt CB -b ${outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv
rm -r ${outs}/${library}/.snakemake
echo ""
echo "Starting AMULET doublet detection"
echo ""
apptainer run -B /fast ${container} AMULET ${outs}/${library}/outs/fragments.tsv.gz ${outs}/${library}/outs/singlecell.csv /opt/AMULET/human_autosomes.txt /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed ${outs}/${library}/AMULET /opt/AMULET/
echo ""
echo "Starting donor SNP genotyping"
echo ""
apptainer run -B /fast ${container} cellsnp-lite -s ${outs}/${library}/outs/possorted_bam.bam -b ${outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv -O ${outs}/${library}/vireo -R /opt/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores --UMItag None
echo ""
echo "Demultiplexing donors with vireo"
echo ""
apptainer run -B /fast ${container} vireo -c ${outs}/${library}/vireo -o ${outs}/${library}/vireo -N $n_donors -p \$num_cores
EOF
            fi
        elif [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            read -p "Would you like to genotype ${library}? (Y/N): " perform_function
            # Convert input to uppercase for case-insensitive comparison
            perform_function=$(echo "$perform_function" | tr '[:lower:]' '[:upper:]')
            # Check if the input is 'N' or 'n'
            if [ "$perform_function" != "Y" ]; then
                echo "Skipping genotyping"
            else
                mkdir -p ${outs}/${library}/AMULET
                mkdir -p ${outs}/${library}/vireo/
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_QC
#SBATCH --output $outs/logs/${library}_genotyping.out
#SBATCH --error $outs/logs/${library}_genotyping.out
#SBATCH --ntasks=32
#SBATCH --mem=68000
#SBATCH --time=8:00:00
num_cores=\$(nproc)
cd ${outs}/${library}
echo "Starting mgatk mtDNA genotyping"
echo ""
apptainer run -B /fast,/usr ${container} mgatk tenx -i ${outs}/${library}/outs/possorted_bam.bam -n output -o ${outs}/${library}/mgatk -c 8 -bt CB -b ${outs}/${library}/outs/filtered_peak_bc_matrix/barcodes.tsv
rm -r ${outs}/${library}/.snakemake
echo ""
echo "Starting AMULET doublet detection"
echo ""
apptainer run -B /fast ${container} AMULET ${outs}/${library}/outs/fragments.tsv.gz ${outs}/${library}/outs/singlecell.csv /opt/AMULET/human_autosomes.txt /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed ${outs}/${library}/AMULET /opt/AMULET/
echo ""
EOF
            fi
        fi
        if [[ "${library}" == *ASAP* ]]; then
            echo "Library ${library} is an ASAP run, performing ADT counting"
            library_csv=${output_project_dir}/${output_project_id}_scripts/libraries/${library}_ADT.csv
            if [[ -f "${library}_csv" ]]; then
                rm ${library}_csv
            fi
            # Initialize associative array for unique entries
            declare -A unique_entries
            # Loop through project directories
            for project_id in "${project_ids[@]}"; do
                project_fastqs="${prefix}/${project_id}/${project_id}_fastq"
                # Loop through matching directories
                for folder in "${project_fastqs}"/*/outs; do
                    # Extract the modality from the folder name
                    if [[ "$folder" == *"ASAP_DI_ADT"* ]]; then
                        modality="ADT"
                    elif [[ "$folder" == *"ASAP_DI_HTO"* ]]; then
                        modality="HTO"
                    else
                        continue
                    fi

                    # Search for matching files and remove duplicates
                    matching_files=($(find "$folder" -type f -name "${library}_*" | grep -E "${library}_(ADT|HTO).*\.fastq\.gz"))
                    # Check if matching files are found
                    if [ ${#matching_files[@]} -gt 0 ]; then
                        # Process each matching file
                        for file_path in "${matching_files[@]}"; do
                            # Extract the fastq_name by removing suffix until $modality
                            fastq_name=${library}_${modality}
                            # Store the directory name
                            directory=$(dirname "$file_path")
                            # Construct the unique identifier
                            identifier="${directory},${fastq_name}"
                            # Append the identifier to the list of unique entries
                            if [ ! "${unique_entries["$identifier"]}" ]; then
                                unique_entries["$identifier"]=1
                            fi
                        done
                    fi
                done
            done
            # Print the unique entries to CSV file
            for entry in "${!unique_entries[@]}"; do
                echo "$entry" >> "${library}_csv"
            done
            # Read the input CSV file line by line
            while IFS= read -r line; do
                # Split the line by comma
                IFS=',' read -r -a parts <<< "$line"
                # Extract fastq_dir and fastq_library
                fastq_dir="${parts[0]}"
                fastq_library="${parts[1]}"
                # Append to the respective variables
                fastq_dirs="$fastq_dirs,$fastq_dir"
                fastq_librarys="$fastq_librarys,$fastq_library"
            done < "${library}_csv"
            # Remove leading comma
            fastq_dirs="${fastq_dirs:1}"
            fastq_librarys="${fastq_librarys:1}"
                echo ""
                echo "For ${library}, the following ASAP FASTQ files will be converted to KITE-compatible FASTQ files "
                echo $fastq_librarys
                echo "In the directories"
                echo $fastq_dirs
                echo ""
                # Ask the user if they want to submit the indices for FASTQ generation
                echo "Do you want to submit with these options? (Y/N)"
                read -r choice
                # Process choices
                if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
                        :
                elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                        exit 1
                else
                        echo "Invalid choice. Exiting"
                fi
                ATAC_whitelist=$HOME/group/work/bin/whitelists/737K-cratac-v1.txt
                ADT_outs=${outs}/${library}/ADT
                mkdir -p ${ADT_outs}
                ADT_index_folder=${outs}/${library}/ADT_index
                mkdir -p ${ADT_index_folder}/temp
                corrected_fastq=${output_project_dir}/${output_project_id}_fastq/ASAP_DI_ADT_corrected
                mkdir -p $corrected_fastq
                ADT_file=${output_project_dir}/${output_project_id}_scripts/ADT_files/${ADT_file}
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}_ADT
#SBATCH --output $outs/logs/${library}_ADT.out
#SBATCH --error $outs/logs/${library}_ADT.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=12:00:00
num_cores=\$(nproc)
cd ${outs}/${library}
echo "Running featuremap"
echo ""
apptainer run -B /fast ${container} featuremap ${ADT_file} --t2g ${ADT_index_folder}/FeaturesMismatch.t2g --fa ${ADT_index_folder}/FeaturesMismatch.fa --header --quiet
echo "Running kallisto index"
echo ""
apptainer run -B /fast ${container} kallisto index -i ${ADT_index_folder}/FeaturesMismatch.idx -k 15 ${ADT_index_folder}/FeaturesMismatch.fa
echo "Running asap_to_kite"
echo ""
apptainer run -B /fast ${container} ASAP_to_KITE -f $fastq_dirs -s $fastq_librarys -o ${corrected_fastq}/${library} -c \$num_cores
echo "Running kallisto bus"
echo ""
apptainer run -B /fast ${container} kallisto bus -i ${ADT_index_folder}/FeaturesMismatch.idx -o ${ADT_index_folder}/temp -x 0,0,16:0,16,26:1,0,0 -t \$num_cores ${corrected_fastq}/${library}*
echo "Running bustools correct"
echo ""
apptainer run -B /fast ${container} bustools correct -w ${whitelist} ${ADT_index_folder}/temp/output.bus -o ${ADT_index_folder}/temp/output_corrected.bus
echo "Running bustools sort"
echo ""
apptainer run -B /fast ${container} bustools sort -t \$num_cores -o ${ADT_index_folder}/temp/output_sorted.bus ${ADT_index_folder}/temp/output_corrected.bus
echo "Running bustools count"
echo ""
apptainer run -B /fast ${container} bustools count -o ${outs}/${library}/ADT/ --genecounts -g ${ADT_index_folder}/FeaturesMismatch.t2g -e ${ADT_index_folder}/temp/matrix.ec -t ${ADT_index_folder}/temp/transcripts.txt ${ADT_index_folder}/temp/output_sorted.bus
rm -r ${ADT_index_folder}
EOF
        fi
        
    else
        # Action when neither file is found
        echo "Neither feature matrix nor peak matrix was found for ${library}"
        exit 1
    fi
done
