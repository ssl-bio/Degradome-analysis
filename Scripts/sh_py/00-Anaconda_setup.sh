#!/bin/bash

# Creates 4 different anaconda environments: 'pydeg_map' 'pydeg_R' 'pydeg_run' and 'pydeg_aux'. The first environment is used for most of the steps including  mapping of the reads; the second is used mainly for running the R scripts. The third environment is used for those programs that require python 2 mainly, PyDegradome and mirmap. The fourth environment is used to install packages that cause conflicts when installed along others.

# Execution: /bin/bash ./Scripts/sh_py/00-Anaconda_setup.sh
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

conda install -c bioconda bioawk fastqc trim-galore bedtools seqkit bowtie2 picard samtools biopython qualimap htseq deeptools -y
conda install -c bioconda agat -y
conda install -c anaconda pandas pillow tk -y
conda install -c conda-forge dash dash-bootstrap-components pdf2image -y
conda install -c plotly plotly plotly_express -y
pip install dash-dangerously-set-inner-html dash_bootstrap_templates nbib dash_breakpoints

set +u
conda deactivate

# Create R environment
conda create --name pydeg_R python=3.11.4 -y
source activate pydeg_R
set -u

# conda install -c conda-forge --file pydeg_R_conda-forge.txt
conda install -c conda-forge pkg-config=0.29.2 r-curl=5.0.1 r-base=4.3.0 -y

set +u
conda deactivate

# Create environment for running pydegradome
conda create --name pydeg_run python=2 -y
source activate pydeg_run
set -u

conda install -c bioconda biopython tophat viennarna unittest2 dendropy -y
conda install -c anaconda scipy pandas pillow tk -y

set +u
conda deactivate

# Create auxiliary environment
conda create --name pydeg_aux python=3 -y
source activate pydeg_aux
set -u

conda install -c bioconda salmon=1.10.3 -y
