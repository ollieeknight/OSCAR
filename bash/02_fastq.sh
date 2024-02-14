#!/bin/bash

# usage is: bash 02_fastq.sh --project-id 'project_id'
# such as bash 02_fastq.sh --project-id K001 --prefix $HOME/scratch/ngs/

# make sure your processing folder structure is similar to follows:

# |── {project_id}_bcl/
# │   ├── Config/
# │   ├── CopyComplete.txt
# │   ├── Data/
# │   ├── InterOp/
# │   ├── Logs/
# │   ├── Recipe/
# │   ├── RTA3.cfg
# │   ├── RTAComplete.txt
# │   ├── RunInfo.xml
# │   ├── RunParameters.xml
# │   ├── SequenceComplete.txt
# │   └── Thumbnail_Images/
# └── {project_id}_scripts/
#     ├── adt_files/
#     ├── indices/
#     └── metadata/ # METADATA MUST BE IN THIS FOLDER!

# Define default values
OSCAR_script_dir=$(dirname "${BASH_SOURCE[0]}")
prefix="$HOME/scratch/ngs"

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project-id)
      project_id="$2"
      shift 2
      ;;
    --prefix)
      prefix="$2"
      shift 2
      ;;
    *)      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

# Check if project_id is empty
if [ -z "$project_id" ]; then
    echo -e "\033[0;31mERROR:\033[0m Please provide a project_id using the --project-id option."
    exit 1
fi

read_length=$(awk -F '"' '/<Read Number="1"/ {print $4}' ${prefix}/${project_id}/${project_id}_bcl/RunInfo.xml)
if [ "$read_length" -gt 45 ]; then
    run_type="ATAC"
elif [ "$read_length" -lt 45 ]; then
    run_type="GEX"
else
    echo -e "\033[0;31mERROR:\033[0m Cannot determine run type, please check ${project_dir}/${project_id}_bcl/RunInfo.xml"
    exit 1
fi

echo ""
echo -e "\033[34mINFO:\033[0m $project_id is an $run_type run, processing appropriately"

project_dir="${prefix}/$project_id"
script_dir="${project_dir}/${project_id}_scripts"
fastq_dir=${project_dir}/${project_id}_fastq
mkdir -p $fastq_dir/logs

# Check if indices folder does not exist, and exit if it's not present
if [ ! -d "$script_dir/indices" ]; then
    echo -e "\033[0;31mERROR:\033[0m Indices folder does not exist, please run the process_metadata script"
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

# Define base masks
xml_file=$project_dir/${project_id}_bcl/RunInfo.xml
if [ ! -f "${xml_file}" ]; then
  echo "Sequencing run RunInfo.xml file not found under $project_dir/${project_id}_bcl/. This is required to determine base masks"
  exit 1
fi
if [[ -f "$xml_file" ]]; then
    if grep -q '<Reads>' "$xml_file" && grep -q '</Reads>' "$xml_file"; then
        num_reads=$(grep -o '<Read Number="' "$xml_file" | wc -l)
        if [ "$num_reads" -eq 3 ]; then
            reads=3
            elif [ "$num_reads" -eq 4 ]; then
                reads=4
            fi
    fi
else
    echo -e "\033[0;31mERROR:\033[0m RunInfo.xml contains unexpected reads, please check the file"
    exit 1
fi

# Define base reads, essentially for whether it is a single or dual index sequencing run
if [[ $reads == 3 ]]; then
        base_mask_SI_3prime_GEX='Y28n*,I8n*,Y90n*'
        base_mask_DI_3prime_GEX='Y28n*,I8n*,Y90n*'
        base_mask_SI_3prime_ADT='Y28n*,I8n*,Y90n*'
        base_mask_DOGMA_ADT='Y28n*,I8n*,Y90n*'
elif [[ $reads == 4 ]]; then
        base_mask_SI_3prime_GEX='Y28n*,I8n*,N*,Y90n*'
        base_mask_SI_5prime_GEX='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DI_3prime_GEX='Y28n*,I8n*,I8n*,Y90n*'
        base_mask_DI_5prime_GEX='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_SI_3prime_ADT='Y28n*,I8n*,N*,Y90n*'
        base_mask_DI_5prime_ADT='Y26n*,I10n*,I10n*,Y90n*'
        base_mask_DOGMA_GEX='Y28n*,I10n*,I10n*,Y90n*'
        base_mask_DOGMA_ATAC='Y100n*,I8n*,Y24n*,Y100n*'
        base_mask_DOGMA_ADT='Y28n*,I8n*,N*,Y90n*'
        base_mask_ASAP_ATAC='Y100n*,I8n*,Y16n*,Y100n*'
        base_mask_ASAP_ADT='Y100n*,I8n*,Y16n*,Y100n*'
else
    echo -e "\033[0;31mERROR:\033[0m Cannot determine number of reads, check RunInfo.xml file"
        exit 1
fi

# List the index files to demultiplex
index_files=($(ls "$project_dir/${project_id}_scripts/indices"))

# Iterate over the index files
for file in "${index_files[@]}"; do
    echo ""
    echo "-------------"
    echo ""
    # Extract the file name without extension
    index_file="${file%.*}"
    if [[ $file == CITE* ]]; then
        cellranger_command='cellranger mkfastq'
        if [[ $file == *_SI_* ]]; then
            index_type='SI'
            filter_option='--filter-single-index'
            if [[ $file == *_GEX* ]]; then
                base_mask=$base_mask_SI_3prime_GEX
            elif [[ $file == *_ADT* ]] || [[ $file == *_HTO* ]]; then
                base_mask=$base_mask_SI_3prime_ADT
            fi
        elif [[ $file == *_DI_* ]]; then
            index_type='DI'
            filter_option='--filter-dual-index'
            if [[ $file == *_3prime* ]]; then
                if [[ $file == *_GEX* ]]; then
                    base_mask=$base_mask_DI_3prime_GEX
                elif [[ $file == *_ADT* ]] || [[ $file == *_HTO* ]]; then
                    base_mask=$base_mask_DI_3prime_ADT
                fi
            elif [[ $file == *_5prime* ]]; then
                if [[ $file == *_GEX* ]] || [[ $file == *_VDJ-T* ]] || [[ $file == *_VDJ-B* ]]; then
                    base_mask=$base_mask_DI_5prime_GEX
                elif [[ $file == *_ADT* ]] || [[ $file == *_HTO* ]]; then
                    base_mask=$base_mask_DI_5prime_ADT
                fi
            fi
        fi
    elif [[ $file == GEX* ]]; then
        cellranger_command='cellranger mkfastq'
        if [[ $file == *_SI* ]]; then
        index_type='SI'
            filter_option='--filter-single-index'
            if [[ $file == *_GEX* ]]; then
                base_mask=$base_mask_SI_3prime_GEX
            fi
        elif [[ $file == *_DI* ]]; then
        index_type='DI'
        filter_option='--filter-dual-index'
            if [[ $file == *_3prime* ]]; then
                if [[ $file == *_GEX* ]]; then
                    base_mask=$base_mask_DI_3prime_GEX
                fi
            elif [[ $file == *_5prime* ]]; then
                if [[ $file == *_GEX* ]]; then
                    base_mask=$base_mask_DI_5prime_GEX
                fi
            fi
        fi
    elif [[ $file == DOGMA* ]]; then
        index_type='DI'
        filter_option='--filter-dual-index'
        if [[ $file == *_GEX* ]]; then
            cellranger_command='cellranger mkfastq'
            base_mask=$base_mask_DOGMA_GEX
        elif [[ $file == *_ADT* ]] || [[ $file == *_HTO* ]]; then
                base_mask=$base_mask_DOGMA_ADT
                if grep -q "ATAC" <<< "$csv_file"; then
                        cellranger_command='cellranger-atac mkfastq'
                        break
                else
                        cellranger_command='cellranger mkfastq'
                fi
        elif [[ $file == *_ATAC* ]]; then
            cellranger_command='cellranger-atac mkfastq'
            base_mask=$base_mask_DOGMA_ATAC
        fi
    elif [[ $file == ASAP* ]]; then
        cellranger_command='cellranger-atac mkfastq'
        index_type='DI'
        filter_option='--filter-dual-index'
        if [[ $file == *_ATAC* ]]; then
            base_mask=$base_mask_ASAP_ATAC
        elif [[ $file == *_ADT* ]] || [[ $file == *_HTO* ]]; then
            base_mask=$base_mask_ASAP_ADT
        fi
    else
        echo -e "\033[0;31mERROR:\033[0m Cannot determine base mask for $index_file, please check path"
        exit 1
    fi

    echo "For index file ${index_file}, the following FASTQ demultiplexing script will be run:"
    echo ""
    echo "$cellranger_command --id ${index_file} --run $project_dir/${project_id}_bcl --csv $script_dir/indices/$file --use-bases-mask $base_mask --delete-undetermined --barcode-mismatches 1 $filter_option"
    echo ""

    # Ask the user if they want to submit the indices for FASTQ generation
    echo -e "\033[0;33mINPUT REQUIRED:\033[0m Is this alright? (Y/N)"
    read -r choice

    # Process choices
    if [ "$choice" = "Y" ] || [ "$choice" = "y" ]; then
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${project_id}
#SBATCH --output $fastq_dir/logs/${index_file}_FASTQ.out
#SBATCH --error $fastq_dir/logs/${index_file}_FASTQ.out
#SBATCH --ntasks=32
#SBATCH --mem=64000
#SBATCH --time=3:00:00
cd $fastq_dir
echo ""
echo "$cellranger_command --id ${index_file} --run $project_dir/${project_id}_bcl --csv $script_dir/indices/$file --use-bases-mask $base_mask --delete-undetermined --barcode-mismatches 1 $filter_option"
echo ""
apptainer run -B /fast ${container} $cellranger_command --id ${index_file} --run ${project_dir}/${project_id}_bcl --csv ${script_dir}/indices/${file} --use-bases-mask ${base_mask} --delete-undetermined --barcode-mismatches 1 ${filter_option}
mkdir -p ${fastq_dir}/${index_file}/fastqc
find "${index_file}/outs/fastq_path/H"* -name "*.fastq.gz" | parallel -j 16 "apptainer run -B /fast ${container} fastqc {} --outdir ${index_file}/fastqc"
apptainer run -B /fast ${container} multiqc "${fastq_dir}/${index_file}" -o "${fastq_dir}/${index_file}/multiqc"
rm -r ${fastq_dir}/${index_file}/_* ${fastq_dir}/${index_file}/MAKE*
EOF
    elif [ "$choice" = "N" ] || [ "$choice" = "n" ]; then
        :
    else
        echo -e "\033[0;31mERROR:\033[0m Invalid choice. Exiting"
    fi
done
