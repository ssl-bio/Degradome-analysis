#!/bin/bash

# Creates 3 different anaconda environments: 'pydeg_map' 'pydeg_R' and 'pydeg_run' The first 2 use a python 3 and the last uses python 2. As they names suggest the first is used mainly for mapping whereas, the second is used exclusively for running the R scripts. The last environment is used for those programs that require python 2 namely, PyDegradome and mirmap.

# Execution: ./Scripts/sh_py/00-Anaconda_setup.sh
#==================================================
eval "$(conda shell.bash hook)"

# Add additional channels (Comment out if already selected or they will take priority)
conda config --add channels conda-forge
conda config --add channels bioconda

# Create mapping environment
conda create --name pydeg_map python=3 -y
source activate pydeg_map

#Abort if get any error
set -eo pipefail

conda install -c bioconda bioawk fastqc trim-galore bedtools seqkit bowtie2 picard samtools biopython qualimap htseq deeptools salmon -y
conda install -c bioconda agat -y
conda install -c anaconda pandas pillow tk -y

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

conda install -c bioconda tophat viennarna unittest2 dendropy -y
conda install -c anaconda pandas pillow tk -y

