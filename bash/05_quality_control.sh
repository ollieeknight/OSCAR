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
#     ├── adt_files/
#     ├── indices/
#     ├── libraries/
#     └── metadata/ # METADATA MUST BE IN THIS FOLDER!

# Parse command-line arguments

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --project-id)
        project_ids_in="$2"
        shift 2
        ;;
        *)
        # Unknown option
        echo "Unknown option: $1"
        exit 1
        ;;
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
output_project_dir=$HOME/scratch/ngs/$output_project_id
outs_folder=$output_project_dir/${output_project_id}_outs

# Check if the output_project_id folder exists
if [ ! -d "$outs_folder" ]; then
    echo "Error: Outs folder ($outs_folder) does not exist. Make sure cellranger has run"
    exit 1
fi

OSCAR_QC=$HOME/group/work/bin/OSCAR/OSCAR_QC.sif

# Check if the output_project_id folder exists
if [ ! -f "$OSCAR_QC" ]; then
    echo "Error: OSCAR container not present, please check $HOME/group/work/bin/OSCAR/OSCAR_container.sif exists"
    exit 1
fi

SNP_list=$HOME/group/work/ref/vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
AMULET_scripts=$HOME/group/work/bin/AMULET

samples=($(find ~/scratch/ngs/$output_project_id/${output_project_id}_outs/ -maxdepth 1 -mindepth 1 -type d -not -name 'logs' -exec basename {} \;))

for sample in "${samples[@]}"; do

        echo ""
        echo "----"
        echo ""

    # Extract variables from the sample name
    assay="${sample%%_*}"
    remainder="${sample#*_}"
    experiment_id="${remainder%%_exp*}"
    remainder="${remainder#*_}"
    remainder="${remainder#*exp}"
    historical_number="${remainder%%_lib*}"
    remainder="${remainder#*lib}"
    replicate="${remainder}"

    # Loop through each project_id
    for project_id in "${project_ids[@]}"; do
        metadata_file="$HOME/scratch/ngs/${project_id}/${project_id}_scripts/metadata/metadata.csv"

        # Check if the metadata file exists
        if [ -f "$metadata_file" ]; then
        echo "Searching $project_id for ${sample}"
            # Read metadata from the CSV file line by line
            while IFS=, read -r -a fields; do
                # Check if all individual fields match the criteria
                if [[ "${fields[0]}" == "$assay" && "${fields[1]}" == "$experiment_id" && "${fields[2]}" == "$historical_number" && "${fields[3]}" == "$replicate" ]]; then
                    n_donors="${fields[9]}"
                    adt_file="${fields[10]}"
            break  # Stop searching once a match is found
                fi
            done < "$metadata_file"
        else
            echo "Error: Metadata file not found for project_id: $project_id"
            exit 1
        fi
    done

    feature_matrix_path=$(find "${outs_folder}/$sample/" -type f -name "raw_feature_bc_matrix.h5" -print -quit)
    peak_matrix_path=$(find "${outs_folder}/$sample/" -type f -name "raw_peak_bc_matrix.h5" -print -quit)

    if [ -n "$feature_matrix_path" ]; then
        if [ -d "${outs_folder}/$sample/cellbender/" ]; then
            echo "Cellbender looks to already have been run for $sample"
        else
            echo "Submitting Cellbender QC for $sample"
            mkdir -p ${outs_folder}/$sample/cellbender/
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_cellbender
#SBATCH --output ${outs_folder}/${sample}/QC.log
#SBATCH --error ${outs_folder}/${sample}/QC.log
#SBATCH --ntasks 1
#SBATCH --partition "gpu"
#SBATCH --gres gpu:1
#SBATCH --cpus-per-task 16
#SBATCH --mem 64000
#SBATCH --time 4:00:00
cd ${outs_folder}/$sample
apptainer run --nv -B /fast ${OSCAR_QC} cellbender remove-background --cuda --input ${feature_matrix_path} --output ${outs_folder}/$sample/cellbender/output.h5
rm ckpt.tar.gz
EOF
        )
        job_id=$(echo "$job_id" | awk '{print $4}')
        echo ""
        fi

        if [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            echo "Skipping genotyping for $sample, as this is either a mouse run, or only contains 1 donor"
            job_id=""
        elif [ -d "${outs_folder}/$sample/vireo/" ]; then
                echo "Genotyping is either ongoing or has finished for $sample"
                job_id=""
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" != "" ]]; then
            echo "Submitting vireo genotyping for $sample"
            mkdir -p ${outs_folder}/$sample/vireo/
            sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_vireo
#SBATCH --output ${outs_folder}/${sample}/QC.log
#SBATCH --error ${outs_folder}/${sample}/QC.log
#SBATCH --ntasks=32
#SBATCH --mem=32000
#SBATCH --time=6:00:00
# The following line ensures that this job runs after the previous job with ID $job_id
#SBATCH --dependency=afterok:$job_id
num_cores=\$(nproc)
cd ${outs_folder}/$sample
apptainer run -B /fast ${OSCAR_QC} cellsnp-lite -s ${outs_folder}/${sample}/outs/per_sample_outs/$sample/count/sample_alignments.bam -b ${outs_folder}/${sample}/cellbender/output_cell_barcodes.csv -O ${outs_folder}/${sample}/vireo -R $SNP_list --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores
apptainer run -B /fast ${OSCAR_QC} vireo -c ${outs_folder}/${sample}/vireo -o ${outs_folder}/${sample}/vireo -N $n_donors -p \$num_cores
EOF
            job_id=""
        elif [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' && "$job_id" == "" ]]; then
            echo "Submitting vireo genotyping for $sample"
            mkdir -p ${outs_folder}/$sample/vireo/
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_vireo
#SBATCH --output ${outs_folder}/${sample}/QC.log
#SBATCH --error ${outs_folder}/${sample}/QC.log
#SBATCH --ntasks=32
#SBATCH --mem=32000
#SBATCH --time=6:00:00
num_cores=\$(nproc)
cd ${outs_folder}/$sample
apptainer run -B /fast ${OSCAR_QC} cellsnp-lite -s ${outs_folder}/${sample}/outs/per_sample_outs/$sample/count/sample_alignments.bam -b ${outs_folder}/${sample}/cellbender/output_cell_barcodes.csv -O ${outs_folder}/${sample}/vireo -R $SNP_list --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores
apptainer run -B /fast ${OSCAR_QC} vireo -c ${outs_folder}/${sample}/vireo -o ${outs_folder}/${sample}/vireo -N $n_donors -p \$num_cores
EOF
            job_id=""
        else
            echo "Cannot determine the number of donors for $sample"
            exit 1
        fi

    elif [ -n "$peak_matrix_path" ]; then

        mkdir -p ${outs_folder}/$sample/amulet
        if [[ "$n_donors" != '0' && "$n_donors" != '1' && "$n_donors" != 'NA' ]]; then
            echo "Submitting ATAC QC for $sample, genotyping $n_donors donors"
            mkdir -p ${outs_folder}/$sample/vireo/
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_QC
#SBATCH --output ${outs_folder}/${sample}/QC.log
#SBATCH --error ${outs_folder}/${sample}/QC.log
#SBATCH --ntasks=32
#SBATCH --mem=66000
#SBATCH --time=8:00:00
num_cores=\$(nproc)
cd ${outs_folder}/$sample
echo "Starting mgatk mtDNA genotyping"
echo ""
apptainer run -B /fast,/usr ${OSCAR_QC} mgatk tenx -i ${outs_folder}/$sample/outs/possorted_bam.bam -n output -o ${outs_folder}/$sample/mgatk -c 8 -bt CB -b ${outs_folder}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv
rm -r ${outs_folder}/$sample/.snakemake
echo ""
echo "Starting AMULET doublet detection"
echo ""
apptainer run -B /fast ${OSCAR_QC} AMULET ${outs_folder}/$sample/outs/fragments.tsv.gz ${outs_folder}/$sample/outs/singlecell.csv /opt/AMULET/human_autosomes.txt /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed ${outs_folder}/$sample/amulet /opt/AMULET/
echo ""
echo "Starting donor SNP genotyping"
echo ""
apptainer run -B /fast ${OSCAR_QC} cellsnp-lite -s ${outs_folder}/$sample/outs/possorted_bam.bam -b ${outs_folder}/$sample/outs/filtered_peak_bc_matrix/barcodes.tsv -O ${outs_folder}/$sample/vireo -R $SNP_list --minMAF 0.1 --minCOUNT 20 --gzip -p \$num_cores --UMItag None
echo ""
echo "Demultiplexing donors with vireo"
echo ""
apptainer run -B /fast ${OSCAR_QC} vireo -c ${outs_folder}/$sample/vireo -o ${outs_folder}/$sample/vireo -N $n_donors -p \$num_cores
EOF
        elif [[ "$n_donors" == '0' || "$n_donors" == '1' || "$n_donors" == 'NA' ]]; then
            echo "Skipping genotyping for $sample, as this is either a mouse run, or only contains 1 donor"
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${experiment_id}_QC
#SBATCH --output ${outs_folder}/${sample}/QC.log
#SBATCH --error ${outs_folder}/${sample}/QC.log
#SBATCH --ntasks=32
#SBATCH --mem=66000
#SBATCH --time=8:00:00
num_cores=\$(nproc)
cd ${outs_folder}/$sample
echo "Starting mgatk mtDNA genotyping"
echo ""
apptainer run -B /fast,/usr ${OSCAR_QC} mgatk tenx -i ${outs_folder}/$sample/outs/possorted_bam.bam -n output -o ${outs_folder}/$sample/mgatk -c 8 -bt CB -b ${outs_folder}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv
rm -r ${outs_folder}/$sample/.snakemake
echo ""
echo "Starting AMULET doublet detection"
echo ""
apptainer run -B /fast ${OSCAR_QC} AMULET ${outs_folder}/$sample/outs/fragments.tsv.gz ${outs_folder}/$sample/outs/singlecell.csv /opt/AMULET/human_autosomes.txt /opt/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed ${outs_folder}/$sample/amulet /opt/AMULET/
echo ""
EOF

        elif [ -d "${outs_folder}/$sample/vireo/" ]; then
            echo "Genotyping is either ongoing or has finished for $sample"
        fi

        if [[ "$sample" == *ASAP* ]]; then
            echo "Sample $sample is an ASAP run, performing ADT counting"

            library_csv=${output_project_dir}/${output_project_id}_scripts/libraries/${sample}_ADT.csv
		rm $library_csv
            declare -A unique_entries

            # Loop through project directories
            for project_id in "${project_ids[@]}"; do
                project_fastqs="$HOME/scratch/ngs/${project_id}/${project_id}_fastq"

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

                    # Search for matching files
                    matching_files=($(find "$folder" -type f -name "${sample}_*" | grep -E "${sample}_(ADT|HTO).*\.fastq\.gz"))

                    # Check if matching files are found
                    if [ ${#matching_files[@]} -gt 0 ]; then
                        for fastq_file in "${matching_files[@]}"; do
                            directory=$(dirname "$fastq_file")
                            identifier="${directory},${sample}_${modality}"
                            if [ ! "${unique_entries["$identifier"]}" ]; then
                                # Print the identifier to CSV file
                                echo "$identifier" >> "$library_csv"
                                # Mark the entry as encountered
                                unique_entries["$identifier"]=1
                            fi
                        done
                    fi
                done
            done

            # Read the input CSV file line by line
            while IFS= read -r line; do
                # Split the line by comma
                IFS=',' read -r -a parts <<< "$line"

                # Extract fastq_dir and fastq_sample
                fastq_dir="${parts[0]}"
                fastq_sample="${parts[1]}"

                # Append to the respective variables
                fastq_dirs="$fastq_dirs,$fastq_dir"
                fastq_samples="$fastq_samples,$fastq_sample"
            done < "$library_csv"

            # Remove leading comma
            fastq_dirs="${fastq_dirs:1}"
            fastq_samples="${fastq_samples:1}"

            # Print the variables
            echo "Fastq Directories: $fastq_dirs"
            echo "Fastq Samples: $fastq_samples"

            whitelist=$HOME/group/work/bin/whitelists/737K-cratac-v1.txt
            adt_outs=${outs_folder}/${sample}/ADT
            mkdir -p ${adt_outs}
            adt_index_folder=${outs_folder}/${sample}/adt_index
            mkdir -p ${adt_index_folder}/temp
            corrected_fastq=${output_project_dir}/${output_project_id}_fastq/ASAP_DI_ADT_corrected
            mkdir -p $corrected_fastq
            adt_file=${output_project_dir}/${output_project_id}_scripts/adt_files/${adt_file}
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${sample}_ADT
#SBATCH --output ${outs_folder}/${sample}/ADT/count.log
#SBATCH --error ${outs_folder}/${sample}/ADT/count.log
#SBATCH --ntasks=32
#SBATCH --mem=64000
#SBATCH --time=12:00:00
num_cores=\$(nproc)
cd ${outs_folder}/$sample
echo "Running featuremap"
echo ""
apptainer run -B /fast ${OSCAR_QC} featuremap ${adt_file} --t2g ${adt_index_folder}/FeaturesMismatch.t2g --fa ${adt_index_folder}/FeaturesMismatch.fa --header --quiet
echo "Running kallisto index"
echo ""
apptainer run -B /fast ${OSCAR_QC} kallisto index -i ${adt_index_folder}/FeaturesMismatch.idx -k 15 ${adt_index_folder}/FeaturesMismatch.fa
echo "Running asap_to_kite"
echo ""
apptainer run -B /fast ${OSCAR_QC} ASAP_to_KITE -f $fastq_dirs -s $fastq_samples -o ${corrected_fastq}/${sample} -c \$num_cores
echo "Running kallisto bus"
echo ""
apptainer run -B /fast ${OSCAR_QC} kallisto bus -i ${adt_index_folder}/FeaturesMismatch.idx -o ${adt_index_folder}/temp -x 0,0,16:0,16,26:1,0,0 -t \$num_cores ${corrected_fastq}/${sample}*
echo "Running bustools correct"
echo ""
apptainer run -B /fast ${OSCAR_QC} bustools correct -w ${whitelist} ${adt_index_folder}/temp/output.bus -o ${adt_index_folder}/temp/output_corrected.bus
echo "Running bustools sort"
echo ""
apptainer run -B /fast ${OSCAR_QC} bustools sort -t \$num_cores -o ${adt_index_folder}/temp/output_sorted.bus ${adt_index_folder}/temp/output_corrected.bus
echo "Running bustools count"
echo ""
apptainer run -B /fast ${OSCAR_QC} bustools count -o ${outs_folder}/${sample}/ADT/ --genecounts -g ${adt_index_folder}/FeaturesMismatch.t2g -e ${adt_index_folder}/temp/matrix.ec -t ${adt_index_folder}/temp/transcripts.txt ${adt_index_folder}/temp/output_sorted.bus
rm -r ${adt_index_folder}
EOF
        fi

    else
        # Action when neither file is found
        echo "Neither feature matrix nor peak matrix was found for $sample"
        exit 1
    fi
done
