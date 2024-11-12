#!/bin/bash

# Enable debugging
set -x

# Default values
oscar_dir=$(dirname "${BASH_SOURCE[0]}")
source "${oscar_dir}/functions.sh"
dir_prefix="${HOME}/scratch/ngs"
metadata_file_name="metadata.csv"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    if [[ "$1" == --* ]]; then
        var_name="${1#--}"
        var_name="${var_name//-/_}"
        declare "$var_name"="$2"
        shift 2
    else
        echo "Invalid option: $1"
        exit 1
    fi
done

check_project_id

# Define project directories
project_dir="${dir_prefix}/${project_id}"
project_scripts="${project_dir}/${project_id}_scripts"
project_indices="${project_scripts}/indices"
project_libraries="${project_scripts}/libraries"
output_project_dir="${project_dir}/output"
fastq_dir=${output_project_dir}/${project_id}_fastq

# Check if metadata file exists
metadata_file="${project_scripts}/metadata/${metadata_file_name}"
check_metadata_file "${metadata_file}"

# Check if indices folder exists
check_folder_exists "${project_scripts}/indices"

# Pull necessary OSCAR containers
check_and_pull_oscar_containers

# List index files and extract flowcell ID from RunInfo.xml
index_files=($(ls "${project_dir}/${project_id}_scripts/indices"))
flowcell_id=$(grep "<Flowcell>" "${project_dir}/${project_id}_bcl/RunInfo.xml" | sed -e 's|.*<Flowcell>\(.*\)</Flowcell>.*|\1|')

# Validate the mode
validate_mode "${mode}"

metadata_file="${project_dir}/${project_id}_scripts/metadata/metadata.csv"

# Take the csv files into a list and remove the .csv suffix
libraries=($(ls "${project_libraries}" | awk -F/ '{print $NF}' | awk -F. '{print $1}'))
outs=${project_dir}/${project_id}_outs/
mkdir -p ${outs}/

# Iterate over each library file to submit counting jobs
for library in "${libraries[@]}"; do
    echo "Processing library ${library}"

if grep -q '.*ADT.*ASAP_.*' "${project_libraries}/${library}.csv"; then
    echo "Library ${library} is an ADT file for ASAP, processing later"
elif grep -q '.*HTO.*ASAP_.*' "${project_libraries}/${library}.csv"; then
    echo "Library ${library} is an HTO file for ASAP, processing later"
    elif grep -q '.*ATAC.*' "${project_libraries}/${library}.csv"; then
    echo "Processing ${library} as an ATAC run"
fastq_names=""
fastq_dirs=""

read fastq_names fastq_dirs < <(count_read_csv "${project_libraries}" "$library")

extra_arguments=$(count_check_dogma "${project_libraries}" "$library")


        echo ""
        echo "cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastq_dirs --sample $fastq_names --localcores $num_cores $extra_arguments"
        echo "(number of cores will change upon submission)"

        # Ask the user if they want to submit the indices for FASTQ generation
        echo -e "\033[0;33mINPUT REQUIRED:\033[0m Is this alright? (Y/N)"
        read -r choice
        while [[ ! $choice =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -r choice
        done
        # Process choices
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
                count_submitted='YES'
            mkdir -p $outs/logs
            # Submit the job to slurm for counting
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/cellranger_${library}.out
#SBATCH --error $outs/logs/cellranger_${library}.out
#SBATCH --ntasks=32
#SBATCH --mem=64000
#SBATCH --time=72:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /data $container cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastq_dirs --sample $fastq_names --localcores \$num_cores $extra_arguments
rm -r $outs/$library/_* $outs/$library/SC_ATAC_COUNTER_CS
EOF
            )
            job_id=$(echo "$job_id" | awk '{print $4}')
        elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
            count_submitted='NO'
        else
            echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
        fi
    if [[ "$library" == *"ASAP"* ]]; then
    # Ask the user if they want to queue ADT counting
    echo "As this is an ASAP-seq run, would you like to queue ADT counting? (Y/N)"
    read -r choice
    while [[ ! $choice =~ ^[YyNn]$ ]]; do
        echo "Invalid input. Please enter Y or N."
        read -r choice
    done
    
    if [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
        continue
    elif [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then

ADT_file=$(count_read_metadata "$metadata_file" "$library")

read fastq_dirs fastq_libraries < <(count_read_adt_csv "${project_libraries}" "$library")

        echo ""
        echo "For ${library}, the following ASAP FASTQ files will be converted to KITE-compatible FASTQ files:"
        echo $fastq_libraries
        echo "In the directories:"
        echo $fastq_dirs
        echo ""
        echo "With the ADT file:"
        echo ${ADT_file}
        echo ""

        # Ask the user if they want to submit with or without dependency
        echo "Do you want to submit with dependency on previous job? (Y/N)"
        read -r choice
        while [[ ! $choice =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -r choice
        done
        
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
            sbatch_dependency="--dependency=afterok:$job_id"
        else
            sbatch_dependency=""
        fi
        
        # Submit the job
        sbatch $sbatch_dependency <<EOF
#!/bin/bash
#SBATCH --job-name ${library}_ADT
#SBATCH --output $outs/logs/kite_${library}.out
#SBATCH --error $outs/logs/kite_${library}.out
#SBATCH --ntasks=16
#SBATCH --mem=96000
#SBATCH --time=72:00:00
num_cores=\$(nproc)
cd ${outs}/${library}
echo "Running featuremap"
echo ""
mkdir -p ${ADT_index_folder}/temp
apptainer run -B /data ${container} featuremap ${ADT_file} --t2g ${ADT_index_folder}/FeaturesMismatch.t2g --fa ${ADT_index_folder}/FeaturesMismatch.fa --header --quiet
echo ""
echo "Running kallisto index"
echo ""
apptainer run -B /data ${container} kallisto index -i ${ADT_index_folder}/FeaturesMismatch.idx -k 15 ${ADT_index_folder}/FeaturesMismatch.fa
echo ""
echo "Running asap_to_kite"
echo ""
mkdir -p $corrected_fastq
apptainer run -B /data ${container} ASAP_to_KITE -f $fastq_dirs -s $fastq_libraries -o ${corrected_fastq}/${library} -c \$num_cores
echo ""
echo "Running kallisto bus"
echo ""
apptainer run -B /data ${container} kallisto bus -i ${ADT_index_folder}/FeaturesMismatch.idx -o ${ADT_index_folder}/temp -x 0,0,16:0,16,26:1,0,0 -t \$num_cores ${corrected_fastq}/${library}*
echo ""
echo "Running bustools correct"
echo ""
apptainer run -B /data ${container} bustools correct -w ${ATAC_whitelist} ${ADT_index_folder}/temp/output.bus -o ${ADT_index_folder}/temp/output_corrected.bus
echo ""
echo "Running bustools sort"
echo ""
apptainer run -B /data ${container} bustools sort -t \$num_cores -o ${ADT_index_folder}/temp/output_sorted.bus ${ADT_index_folder}/temp/output_corrected.bus
echo ""
echo "Running bustools count"
echo ""
mkdir -p ${ADT_outs}
apptainer run -B /data ${container} bustools count -o ${outs}/${library}/ADT/ --genecounts -g ${ADT_index_folder}/FeaturesMismatch.t2g -e ${ADT_index_folder}/temp/matrix.ec -t ${ADT_index_folder}/temp/transcripts.txt ${ADT_index_folder}/temp/output_sorted.bus
rm -r ${ADT_index_folder}
EOF

        fi
        
        # Reset the fastqs variable
        job_id=""
        fastqs=""
    fi
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX*' "${project_libraries}/${library}.csv"; then
    echo "Processing ${library} as an CITE/GEX run"
        echo ""
        echo "For library $library"
        echo ""
        cat ${project_libraries}/${library}.csv
        echo ""
        echo "cellranger multi --id $library --csv ${project_libraries}/${library}.csv --localcores $num_cores"
        echo "(number of cores will change upon submission)"

        # Ask the user if they want to submit the indices for FASTQ generation
        echo -e "\033[0;33mINPUT REQUIRED:\033[0m Is this alright? (Y/N)"
        read -r choice
        while [[ ! $choice =~ ^[YyNn]$ ]]; do
            echo "Invalid input. Please enter Y or N."
            read -r choice
        done
        # Process choices
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
            mkdir -p $outs/logs/
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_counting.out
#SBATCH --error $outs/logs/${library}_counting.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=96:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /data "$container" cellranger multi --id "${library}" --csv "${project_libraries}/${library}.csv" --localcores "\$num_cores"
rm -r $outs/$library/SC_MULTI_CS $outs/$library/_*
EOF
        fi
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine whether ${library} is a GEX or ATAC run. is 'GEX' or 'ATAC present in its library csv file?"
        exit 1
    fi
    echo ""
    echo "-------------"
    echo ""
done

# Disable debugging
set +x