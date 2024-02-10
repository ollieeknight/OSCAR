#!/bin/bash

# usage is: bash 02_fastq.sh --project-id 'project_id'

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
#     ├── metadata/ # METADATA MUST BE IN THIS FOLDER!
#     └── process_metadata.sh

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
    echo "Please provide a project_id using the --project-id option."
    exit 1
fi

project_dir=$HOME/scratch/ngs/$project_id

script_dir=${project_dir}/${project_id}_scripts
fastq_dir=${project_dir}/${project_id}_fastq
mkdir -p $fastq_dir/logs

# Check if indices folder does not exist, and exit if it's not present
if [ ! -d "$script_dir/indices" ]; then
    echo "Indices folder does not exist, please run the process_metadata script"
    exit 1
fi

# Check the symbolic link for the group folder in the users $HOME
if [ ! -d "$HOME/group" ]; then
    ln -sr /fast/work/groups/ag_romagnani/ $HOME/group
fi

container=$HOME/group/work/bin/OSCAR/OSCAR_counting.sif

# Check that the singularity container is available
if [ ! -f "${container}" ]; then
  echo "OSCAR_counting.sif singularity file not found, please check path"
  exit 1
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
        echo "RunInfo.xml contains unexpected reads, please check the file"
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
        echo "Cannot determine number of reads, check RunInfo.xml file"
        exit 1
fi

# List the index files to demultiplex
index_files=($(ls "$project_dir/${project_id}_scripts/indices"))

# Iterate over the index files
for file in "${index_files[@]}"; do
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
            filter_option='--filter-dual-index'
            index_type='DI'
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
            cellranger_command='cellranger mkfastq'
            base_mask=$base_mask_DOGMA_ADT
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
        echo "Cannot determined base mask for $index_file"
        exit 1
    fi

    echo "$cellranger_command --id ${index_file} --run $project_dir/${project_id}_bcl --csv $script_dir/indices/$file --use-bases-mask $base_mask --delete-undetermined --barcode-mismatches 1 $filter_option"
    echo ""
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${index_file}
#SBATCH --output $fastq_dir/logs/${index_file}_SLURM.out
#SBATCH --error $fastq_dir/logs/${index_file}_SLURM.out
#SBATCH --ntasks=32
#SBATCH --mem=64000
#SBATCH --time=3:00:00
cd $fastq_dir
apptainer exec -B /fast,/data ${container} $cellranger_command --id ${index_file} --run ${project_dir}/${project_id}_bcl --csv ${script_dir}/indices/${file} --use-bases-mask ${base_mask} --delete-undetermined --barcode-mismatches 1 ${filter_option}
mkdir -p ${fastq_dir}/${index_file}/fastqc
find "${index_file}/outs/fastq_path/H"* -name "*.fastq.gz" | parallel -j 16 "fastqc {} --outdir ${index_file}/fastqc"
multiqc "${fastq_dir}/${index_file}" -o "${fastq_dir}/${index_file}/multiqc"
rm -r ${fastq_dir}/${index_file}/_* ${fastq_dir}/${index_file}/MAKE*
EOF
    echo ""
    echo "-------------"
    echo ""
done
