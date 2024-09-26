#!/bin/bash

# Define default values
OSCAR_script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
OSCAR_base_dir=$(dirname "$OSCAR_script_dir")
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

metadata_file="${project_dir}/${output_project_id}_scripts/metadata/metadata.csv"

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
    echo "Processing library ${library}"

if grep -q '.*ADT.*ASAP_.*' "${library_folder}/${library}.csv"; then
    echo "Library ${library} is an ADT file for ASAP, processing later"
elif grep -q '.*HTO.*ASAP_.*' "${library_folder}/${library}.csv"; then
    echo "Library ${library} is an HTO file for ASAP, processing later"
    elif grep -q '.*ATAC.*' "${library_folder}/${library}.csv"; then
    echo "Processing ${library} as an ATAC run"
fastq_names=""
fastq_dirs=""

# Read the CSV file line by line
while IFS= read -r line; do
    # Check if the line contains "ATAC"
    if [[ $line == *ATAC* ]]; then
        # Extract fastq name and directory from the line
        IFS=',' read -r fastq_name fastq_dir <<< "$line"

        # Concatenate the fastq name to fastq_names variable
        if [ -n "$fastq_names" ]; then
            fastq_names="${fastq_names},${fastq_name}"
        else
            fastq_names="$fastq_name"
        fi

        # Concatenate the fastq directory to fastq_dirs variable
        if [ -n "$fastq_dirs" ]; then
            fastq_dirs="${fastq_dirs},${fastq_dir}"
        else
            fastq_dirs="$fastq_dir"
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
#SBATCH --output $outs/logs/${library}_counting.out
#SBATCH --error $outs/logs/${library}_counting.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=24:00:00
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
        # Read the metadata file line by line
        while IFS=',' read -r -a fields; do
            assay="${fields[0]}"
            experiment_id="${fields[1]}"
            historical_number="${fields[2]}"
            replicate="${fields[3]}"

            expected_library="${assay}_${experiment_id}_exp${historical_number}_lib${replicate}"

            if [ "$expected_library" == "$library" ]; then
                ADT_file="${fields[10]}"
                break
            fi
        done < "$metadata_file"

        library_csv=${library_folder}/${library}_ADT.csv
        fastq_dirs=''
        fastq_libraries=''
        while IFS= read -r line; do
            IFS=',' read -r -a parts <<< "$line"
            fastq_library="${parts[0]}"
            fastq_dir="${parts[1]}"
            fastq_libraries="$fastq_libraries,$fastq_library"
            fastq_dirs="$fastq_dirs,$fastq_dir"
        done < "${library_csv}"
        fastq_dirs="${fastq_dirs:1}"
        fastq_libraries="${fastq_libraries:1}"
        ATAC_whitelist=${OSCAR_base_dir}/whitelists/737K-cratac-v1.txt
        ADT_outs=${outs}/${library}/ADT
        ADT_index_folder=${outs}/${library}/ADT_index
        corrected_fastq=${project_dir}/${output_project_id}_fastq/ASAP_SI_ADT_corrected
        ADT_file=${project_dir}/${output_project_id}_scripts/ADT_files/${ADT_file}.csv

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
#SBATCH --output $outs/logs/${library}_ADT.out
#SBATCH --error $outs/logs/${library}_ADT.out
#SBATCH --ntasks=16
#SBATCH --mem=128000
#SBATCH --time=12:00:00
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
    elif grep -q '.*GEX*' "${library_folder}/${library}.csv"; then
    echo "Processing ${library} as an CITE/GEX run"
        echo ""
        echo "For library $library"
        echo ""
        cat ${library_folder}/${library}.csv
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
#SBATCH --ntasks=64
#SBATCH --mem=96000
#SBATCH --time=96:00:00
num_cores=\$(nproc)
cd $outs
apptainer run -B /fast,/data "$container" cellranger multi --id "${library}" --csv "${library_folder}/${library}.csv" --localcores "\$num_cores"
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
