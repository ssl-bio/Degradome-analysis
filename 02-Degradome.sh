#!/bin/bash

# Description: Wrapper to further analyze the classification results from PyDegradome (PyDegradome, Gaglia et al. PLOS Pathogens 2015:11) in R.
# 00-Initialization.R: Generates variables and objects used in the rest of the scripts
# 01-Annotation-shared_peak_classification.R: Annotates peaks, calculates signal:noise and classifies shared peaks
# 02-Filtering_Classification_Pooling.R: 
# 03-Drawing_Dplots.R: 
# 04-GetPeakSeq.R : 
# 05-Peak_miRNA_alignment.py: 
# 05-Report-settings.R: 

# Execution: 02-Degradome.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash 02-Degradome.sh Zhang-2021 Zhang-2021_vars.txt

# Check the number of arguments
if [[ -z $1 ]]
then
    echo "No basename for input/output was provided. Exiting..."
    exit 1
else
    echo "$1 will be used as basename for input/output"
    ibase=$1
fi

# Import variables
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/Aux/00-Variable_setup.sh ${1} ${2}
else
    source ${ivars}
fi

#activate conda environment
eval "$(conda shell.bash hook)"
source activate ${conda_pydeg_R}

#Abort if get any error
set -eo pipefail

dir=$(pwd)
# echo $dir

# Run RScripts
echo "Running 00-Initialization.R"
outfile=${supp_data_dir}/"R/Initialization_variables.RData"
if [[ ! -f $outfile ]]
then
    Rscript ${script_dir}/R/00-Initialization.R -d $dir -b $ibase
fi

echo "Running 01-Annotation-shared_peak_classification.R"
Rscript ${script_dir}/R/01-Annotation-shared_peak_classification.R -d $dir -b $ibase

echo "Running 02-Filtering_Classification_Pooling.R"
Rscript ${script_dir}/R/02-Filtering_Classification_Pooling.R -d $dir -b $ibase

echo "Running 03-Drawing_Dplots.R"
Rscript ${script_dir}/R/03-Drawing_Dplots.R -d $dir -b $ibase

echo "Running 04-GetPeakSeq.R "
Rscript ${script_dir}/R/04-GetPeakSeq.R -d $dir -b $ibase

echo "Running 05-Peak_miRNA_alignment.py"
set +u
conda deactivate
source activate ${conda_pydeg_map}
set -u
inputfile=${supp_data_dir}/miRNA_seq/miRNA_sequences.fa
if [[ ! -f ${inputfile} ]]
then
    /bin/bash ${script_dir}/Aux/05-Get_miRNA_seqs.sh $dir $ibase
fi
python ${script_dir}/Aux/05-Peak_miRNA_alignment.py -rd $dir -bn $ibase


echo "Running 05-Report-settings.R"
conda deactivate
source activate ${conda_pydeg_R}
set -u
Rscript ${script_dir}/R/05-Report_main.R -d $dir -b $ibase

