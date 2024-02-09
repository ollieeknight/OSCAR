#!/bin/bash/

cd $PWD
mkdir -p test_run/test_run_bcl
wget -O test_run/cellranger-tiny-bcl-1.2.0.tar.gz https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
tar -xf test_run/cellranger-tiny-bcl-1.2.0.tar.gz -C test_run/test_run_bcl
rm test_run/cellranger-tiny-bcl-1.2.0.tar.gz

mkdir -p test_run/test_run_scripts/indices/
wget -O test_run/test_run_scripts/indices/test_run.csv https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
