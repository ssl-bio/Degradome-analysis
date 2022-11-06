#!/bin/bash

# Description: Wrapper to run R scripts

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
    /bin/bash Scripts/Aux/00-Variable_setup.sh ${1} Degradome_vars_mint.txt
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
Rscript ${script_dir}/R/01-Annotation-shared_peak_classification.R -d $dir

echo "Running 02-Filtering_Classification_Pooling.R"
Rscript ${script_dir}/R/02-Filtering_Classification_Pooling.R -d $dir

echo "Running 03-Drawing_Dplots.R"
Rscript ${script_dir}/R/03-Drawing_Dplots.R -d $dir

echo "Running 04-GetPeakSeq.R "
Rscript ${script_dir}/R/04-GetPeakSeq.R -d $dir

echo "Running 05-Peak_miRNA_alignment.py"
set +u
conda deactivate
source activate ${conda_pydeg_map}
set -u
inputfile=${supp_data_dir}/miRNA_seq/miRNA_sequences.fa
if [[ ! -f ${inputfile} ]]
then
    /bin/bash ${script_dir}/Aux/05-Get_miRNA_seqs.sh $dir
fi
python ${script_dir}/Aux/05-Peak_miRNA_alignment.py -rd $dir


echo "Running 05-Report-settings.R"
conda deactivate
source activate ${conda_pydeg_R}
set -u
Rscript ${script_dir}/R/05-Report-settings.R -d $dir

