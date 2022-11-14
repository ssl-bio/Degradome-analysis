#!/bin/bash

# Script to setup anaconda environment
#==================================================
eval "$(conda shell.bash hook)"

# Add additional channels
conda config --add channels conda-forge
conda config --add channels bioconda

# Create mapping environment
conda create --name pydeg_map python=3 -y
source activate pydeg_map

#Abort if get any error
set -eo pipefail

conda install -c bioconda fastqc trim-galore bowtie2 picard samtools biopython qualimap htseq deeptools salmon seqkit bedtools -y

set +u
conda deactivate

# Create R environment
conda create --name pydeg_R python=3 -y
source activate pydeg_R
set -u

conda install -c conda-forge pkg-config r-curl r-base -y

set +u
conda deactivate
# Create auxiliary environment
conda create --name pydeg_run python=2 -y
source activate pydeg_run

conda install -c bioconda tophat -y
