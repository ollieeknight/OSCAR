#!/bin/bash

# usage is: bash 04_count.sh -project-id project_id

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
        output_project_id="$2"
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
    echo "Please provide a project_id using the --project-id option."
    exit 1
fi

container=$HOME/group/work/bin/OSCAR/OSCAR_counting.sif

if [ ! -f "$container" ]; then
  echo "OSCAR_counting.sig container file not found, please check path"
  exit 1
fi

# Define the project directory
outs="$HOME/scratch/ngs/${output_project_id}/${output_project_id}_outs"
mkdir -p "$outs"

library_folder=$HOME/scratch/ngs/${output_project_id}/${output_project_id}_scripts/libraries

if [ ! -d "$library_folder" ]; then
  echo "Libraries folder not found - did you run the last script (process_libraries.sh)?"
  exit 1
fi

# Take the csv files into a list and remove the .csv suffix
libraries=($(ls "$library_folder" | awk -F/ '{print $NF}' | awk -F. '{print $1}'))

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
        if [ -d "$library/outs" ]; then
            echo "Cellranger looks to have completed for $library, skipping. Check if you do not believe this is the case"
        else
            echo "cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastqs --sample $fastq_name --localcores $num_cores $extra_arguments"
            mkdir -p $outs/logs
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_SLURM.out
#SBATCH --error $outs/logs/${library}_SLURM.out
#SBATCH --ntasks=64
#SBATCH --mem=96000
#SBATCH --time=24:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /fast $container cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $fastqs --sample $fastq --localcores \$num_cores $extra_arguments
rm -r $outs/$library/_* $outs/$library/SC_ATAC_COUNTER_CS
EOF
        # Reset the fastqs variable
        fastqs=""
        fi
    # Check if the modality GEX appears anywhere in the csv file. cellranger multi will process this
    elif grep -q '.*GEX.*' "${library_folder}/${library}.csv"; then
        if [ -d "$library/outs" ]; then
            echo "Cellranger looks to have completed for $library, skipping. Check if you do not believe this is the case"
        else
            echo "cellranger multi --id $library --csv ${library_folder}/${library}.csv"
            mkdir -p $outs/logs/
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_SLURM.out
#SBATCH --error $outs/logs/${library}_SLURM.out
#SBATCH --ntasks=64
#SBATCH --mem=96000
#SBATCH --time=24:00:00
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
