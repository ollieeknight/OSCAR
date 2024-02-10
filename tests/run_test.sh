#!/bin/bash/

# Parse command line arguments using getopts_long function
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --reference)
      reference="$2"
      shift 2
      ;;
    *)      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

# Check if project_id is empty
if [ -z "$reference" ]; then
    echo "Please provide a human reference transcriptome with the --reference option to run this test"
    exit 1
fi


num_cores=$(nproc)

cd $PWD

echo ""
echo "Downloading required run files"
echo ""
mkdir -p test_run/
wget -O test_run/cellranger-tiny-bcl-1.2.0.tar.gz https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz --quiet
tar -xf test_run/cellranger-tiny-bcl-1.2.0.tar.gz -C test_run/
mv test_run/cellranger-tiny-bcl-1.2.0 test_run/test_run_bcl
rm test_run/cellranger-tiny-bcl-1.2.0.tar.gz

mkdir -p test_run/test_run_scripts/indices/
wget -O test_run/test_run_scripts/indices/test_run.csv https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv --quiet

if [ ! -f "oscar-counting_v1.sif" ]; then
	echo "Pulling OSCAR counting image, this might take some time..."
	echo ""
	apptainer pull --arch amd64 library://romagnanilab/default/oscar-counting:v1
fi

cd test_run
echo "Running cellranger mkfastq"
echo ""
apptainer run -B /fast ../oscar-counting_v1.sif cellranger mkfastq --id test_run_fastq --run test_run_bcl/ --csv test_run_scripts/indices/test_run.csv &> mkfastq.log
echo "Running cellranger count"
echo ""
apptainer run -B /fast ../oscar-counting_v1.sif cellranger count --id test_run_sample --fastqs test_run_fastq/outs/fastq_path/H35KCBCXY/ --sample test_sample --localcores $num_cores --transcriptome $reference --chemistry SC3Pv3 --no-bam &> count.log

if [ ! -f "test_run_sample/outs/" ]; then
	echo "Sucess! Shutting down"
else
	echo "Uh oh"
fi
