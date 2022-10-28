#!/bin/bash

# Script to setup anaconda environment
#==================================================
eval "$(conda shell.bash hook)"

# Add additional channels
conda config --add channels conda-forge
conda config --add channels bioconda

# Create main environment
conda create --name pydeg-main python=3 -y
source activate pydeg-main

#Abort if get any error
set -eo pipefail

conda install -c conda-forge pkg-config r-curl r-base -y
conda install -c bioconda fastqc trim-galore bowtie2 picard samtools biopython qualimap htseq deeptools -y


set +u
conda deactivate
source activate $conda_env_py2
set -u

# Create auxiliary environment
conda create --name pydeg-aux python=2 -y
source activate pydeg-aux

conda install -c bioconda tophat -y
