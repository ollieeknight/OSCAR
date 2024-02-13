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
    echo "Please provide a project_id using the --project-id option."
    exit 1
fi

container=$TMPDIR/oscar-counting_latest.sif

# Check that the singularity container is available
if [ ! -f "${container}" ]; then
    echo "oscar-counting_latest.sif singularity file not found, pulling..."
    apptainer pull --dir $TMPDIR library://romagnanilab/default/oscar-counting:latest
fi

# Check that the singularity container is available
if [ ! -f "${container}" ]; then
    echo "oscar-counting_latest.sif singularity file still not found"
        exit 1
fi

# Define the project directory
project_dir=${prefix}/${output_project_id}
library_folder=${project_dir}/${output_project_id}_scripts/libraries
outs="${prefix}/${output_project_id}/${output_project_id}_outs"
mkdir -p "$outs"

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
        fastq_name="${library}_ATAC"
        # Read each line containing "ATAC" in the CSV file, if sequenced across multiple runs
while IFS= read -r line; do
    # Extract the directory part and the fastq_name part
    directory=$(dirname "$line")
    fastq_name=$(basename "$line")

    # If there are multiple directories, merge them with a comma
    if [ "${#directories[@]}" -gt 1 ]; then
        directory=$(IFS=,; echo "${directories[*]}")
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
		echo ""
		echo "For $library, the counting command will be "
		num_cores=$(nproc)
		echo "cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $directory --sample $fastq_name --localcores $num_cores $extra_arguments"
		echo ""
		echo "Where the FASTQ files to input is/are"
		echo "$directory"
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
            mkdir -p $outs/logs
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_counting.out
#SBATCH --error $outs/logs/${library}_counting.out
#SBATCH --ntasks=64
#SBATCH --mem=96000
#SBATCH --time=24:00:00
num_cores=\$(nproc)
cd $outs
echo ""
echo "cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $directory --sample $fastq_name --localcores $num_cores $extra_arguments"
echo ""
apptainer run -B /fast $container cellranger-atac count --id $library --reference $HOME/group/work/ref/hs/GRCh38-hardmasked-optimised-arc/ --fastqs $directory --sample $fastq_name --localcores \$num_cores $extra_arguments
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
                echo ""
                echo "For $library, the counting command will be "
		num_cores=$(nproc)
                echo "cellranger multi --id $library --csv ${library_folder}/${library}.csv --localcores $num_cores"
                echo ""
                echo "Where the library csv input is"
                cat ${library_folder}/${library}.csv
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

            mkdir -p $outs/logs/
            # Submit the job to slurm for counting
            sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${library}
#SBATCH --output $outs/logs/${library}_counting.out
#SBATCH --error $outs/logs/${library}_counting.out
#SBATCH --ntasks=64
#SBATCH --mem=96000
#SBATCH --time=24:00:00
num_cores=\$(nproc)
cd $outs
echo ""
echo "cellranger multi --id $library --csv ${library_folder}/${library}.csv --localcores $num_cores"
echo ""
apptainer run -B /fast $container cellranger multi --id "$library" --csv "${library_folder}/${library}.csv" --localcores \$num_cores
rm -r $outs/$library/SC_MULTI_CS $outs/$library/_*
EOF
        fi
    fi
    echo ""
    echo "-------------"
    echo ""
done
