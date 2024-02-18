#!/bin/bash

# usage is: bash 04_count.sh -project-id project_id
# Requirements at the moment for the reference transcriptomes are
# under here $HOME/group/work/ref/

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

# Define default values
OSCAR_script_dir=$(dirname "${BASH_SOURCE[0]}")
prefix="$HOME/scratch/ngs"
num_cores=$(nproc)

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --project-id)
        output_project_id="$2"
        shift 2
        ;;
        --prefix)
        prefix="$2"
        shift 2
        ;;
        *)
        # Unknown option
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Check if project_id is empty
if [ -z "$output_project_id" ]; then
    echo -e "\033[0;31mERROR:\033[0m Please provide a project_id using the --project-id option."
    exit 1
fi

read_length=$(awk -F '"' '/<Read Number="1"/ {print $4}' ${prefix}/${output_project_id}/${output_project_id}_bcl/RunInfo.xml)
if [ "$read_length" -gt 45 ]; then
    run_type="ATAC"
elif [ "$read_length" -lt 45 ]; then
    run_type="GEX"
else
    echo -e "\033[0;31mERROR:\033[0m Cannot determine run type, please check ${project_dir}/${output_project_id}_bcl/RunInfo.xml"
    exit 1
fi

echo ""
echo -e "\033[34mINFO:\033[0m $output_project_id is an $run_type run, processing appropriately"

project_dir="${prefix}/$output_project_id"
library_folder="${project_dir}/${output_project_id}_scripts/libraries"
fastq_dir=${project_dir}/${output_project_id}_fastq
mkdir -p $fastq_dir/logs

# Check if indices folder does not exist, and exit if it's not present
if [ ! -d "$library_folder" ]; then
    echo -e "\033[0;31mERROR:\033[0m Libraries folder does not exist, please run the process_metadata script"
    exit 1
fi

# Check the symbolic link for the group folder in the users $HOME
if [ ! -d "$HOME/group" ]; then
    ln -s /fast/work/groups/ag_romagnani/ $HOME/group
fi

container=$TMPDIR/oscar-count_latest.sif

# Check that the singularity container is available
if [ ! -f "${container}" ]; then
    echo "oscar-count_latest.sif singularity file not found, pulling..."
    mkdir -p $TMPDIR
    apptainer pull --dir $TMPDIR library://romagnanilab/default/oscar-count:latest
fi

# Take the csv files into a list and remove the .csv suffix
libraries=($(ls "$library_folder" | awk -F/ '{print $NF}' | awk -F. '{print $1}'))
outs=${project_dir}/${output_project_id}_outs/
mkdir -p ${outs}/

# Iterate over each library file to submit counting jobs
for library in "${libraries[@]}"; do
    # If the library file contains the string 'ATAC', it will be counted using cellranger-atac
    if grep -q '.*ATAC.*' "${library_folder}/${library}.csv"; then
        fastq="${library}_ATAC"
        # Read each line containing "ATAC" in the CSV file, if sequenced across multiple runs
        while IFS= read -r line; do
            if [[ $line == *ATAC* ]]; then
                # Concatenate the line with the result, separated by a comma
                if [ -n "$fastq" ]; then
                    fastqs="${fastq},${line}"
                else
                    fastqs="$line"
                fi
            fi
        done < "${library_folder}/${library}.csv"
        # If it's a DOGMA-seq run, chemistry needs to be specified
        if grep -q '.*DOGMA.*' "${library_folder}/${library}.csv"; then
            cat ${library_folder}/${library}.csv
            extra_arguments="--chemistry ARC-v1"
            echo "Adding $extra_arguments as it is a DOGMA/MULTIOME run"
        else
            extra_arguments=""
        fi

        echo ""
        echo "cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastqs --sample $fastq_name --localcores $num_cores $extra_arguments"
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
            mkdir -p $outs/logs
            # Submit the job to slurm for counting
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_counting.out
#SBATCH --error $outs/logs/${library}_counting.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=24:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /fast $container cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastqs --sample $fastq --localcores \$num_cores $extra_arguments
rm -r $outs/$library/_* $outs/$library/SC_ATAC_COUNTER_CS
EOF
            )
            job_id=$(echo "$job_id" | awk '{print $4}')
        elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
            :
        else
            echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
        fi
        if [[ "$library" == *"ASAP"* ]]; then
            # Ask the user if they want to submit the indices for FASTQ generation
            echo "As this is an ASAP-seq run, would you like to queue ADT counting? (Y/N)"
            read -r choice
            while [[ ! $choice =~ ^[YyNn]$ ]]; do
                echo "Invalid input. Please enter Y or N."
                read -r choice
            done
            # Process choices
            if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
                library_csv=${library_folder}/${library}_ADT.csv
                fastq_dirs=''
                fastq_libraries=''
                while IFS= read -r line; do
                    # Split the line by comma
                    IFS=',' read -r -a parts <<< "$line"
                    # Extract fastq_dir and fastq_library
                    fastq_library="${parts[0]}"
                    fastq_dir="${parts[1]}"
                    # Append to the respective variables
                    fastq_libraries="$fastq_libraries,$fastq_library"
                    fastq_dirs="$fastq_dirs,$fastq_dir"
                done < "${library_csv}"
                # Remove leading comma
                fastq_dirs="${fastq_dirs:1}"
                fastq_libraries="${fastq_libraries:1}"
                echo ""
                echo "For ${library}, the following ASAP FASTQ files will be converted to KITE-compatible FASTQ files:"
                echo $fastq_libraries
                echo "In the directories:"
                echo $fastq_dirs
                echo ""
                # Ask the user if they want to submit the indices for FASTQ generation
                echo "Do you want to submit with these options? (Y/N)"
                read -r choice
                while [[ ! $choice =~ ^[YyNn]$ ]]; do
                    echo "Invalid input. Please enter Y or N."
                    read -r choice
                done
                # Process choices
                if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
                    ATAC_whitelist=${OSCAR_script_dir}/../whitelists/737K-cratac-v1.txt
                    ADT_outs=${outs}/${library}/ADT
                    ADT_index_folder=${outs}/${library}/ADT_index
                    corrected_fastq=${project_dir}/${output_project_id}_fastq/ASAP_DI_ADT_corrected
                    ADT_file=${project_dir}/${output_project_id}_scripts/ADT_files/${ADT_file}.csv
sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name ${library}_ADT
#SBATCH --output $outs/logs/${library}_ADT.out
#SBATCH --error $outs/logs/${library}_ADT.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=12:00:00
# The following line ensures that this job runs after the previous job with ID $job_id
#SBATCH --dependency=afterok:$job_id
num_cores=\$(nproc)
cd ${outs}/${library}
echo "Running featuremap"
echo ""
mkdir -p ${ADT_index_folder}/temp
apptainer run -B /fast ${container} featuremap ${ADT_file} --t2g ${ADT_index_folder}/FeaturesMismatch.t2g --fa ${ADT_index_folder}/FeaturesMismatch.fa --header --quiet
echo ""
echo "Running kallisto index"
echo ""
apptainer run -B /fast ${container} kallisto index -i ${ADT_index_folder}/FeaturesMismatch.idx -k 15 ${ADT_index_folder}/FeaturesMismatch.fa
echo ""
echo "Running asap_to_kite"
echo ""
mkdir -p $corrected_fastq
apptainer run -B /fast ${container} ASAP_to_KITE -f $fastq_dirs -s $fastq_libraries -o ${corrected_fastq}/${library} -c \$num_cores
echo ""
echo "Running kallisto bus"
echo ""
apptainer run -B /fast ${container} kallisto bus -i ${ADT_index_folder}/FeaturesMismatch.idx -o ${ADT_index_folder}/temp -x 0,0,16:0,16,26:1,0,0 -t \$num_cores ${corrected_fastq}/${library}*
echo ""
echo "Running bustools correct"
echo ""
apptainer run -B /fast ${container} bustools correct -w ${ATAC_whitelist} ${ADT_index_folder}/temp/output.bus -o ${ADT_index_folder}/temp/output_corrected.bus
echo ""
echo "Running bustools sort"
echo ""
apptainer run -B /fast ${container} bustools sort -t \$num_cores -o ${ADT_index_folder}/temp/output_sorted.bus ${ADT_index_folder}/temp/output_corrected.bus
echo ""
echo "Running bustools count"
echo ""
mkdir -p ${ADT_outs}
apptainer run -B /fast ${container} bustools count -o ${outs}/${library}/ADT/ --genecounts -g ${ADT_index_folder}/FeaturesMismatch.t2g -e ${ADT_index_folder}/temp/matrix.ec -t ${ADT_index_folder}/temp/transcripts.txt ${ADT_index_folder}/temp/output_sorted.bus
rm -r ${ADT_index_folder}
EOF
                elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                    :
                fi
            elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
                :
            fi
        fi
    
        # Reset the fastqs variable
        job_id=""
        fastqs=""
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX.*' "${library_folder}/${library}.csv"; then

        echo ""
        echo "cellranger multi --id $library --csv ${library_folder}/${library}.csv --localcores $num_cores"
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
#SBATCH --time=48:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /fast $container cellranger multi --id "$library" --csv "${library_folder}/${library}.csv" --localcores \$num_cores
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
