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

container=$TMPDIR/oscar-counting_latest.sif

# Check that the singularity container is available
if [ ! -f "${container}" ]; then
    echo "oscar-counting_latest.sif singularity file not found, pulling..."
    mkdir -p $TMPDIR
    apptainer pull --dir $TMPDIR library://romagnanilab/default/oscar-counting:latest
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
        # Process choices
        if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
            mkdir -p $outs/logs
            # Submit the job to slurm for counting
            sbatch <<EOF
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
        elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
            :
        else
            echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
        fi
        # Reset the fastqs variable
        fastqs=""
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX.*' "${library_folder}/${library}.csv"; then

        echo ""
        echo "cellranger multi --id $library --csv ${library_folder}/${library}.csv --localcores $num_cores"
        echo "(number of cores will change upon submission)"

        # Ask the user if they want to submit the indices for FASTQ generation
        echo -e "\033[0;33mINPUT REQUIRED:\033[0m Is this alright? (Y/N)"
        read -r choice
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
    fi
    echo ""
    echo "-------------"
    echo ""
done
