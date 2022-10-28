#!/bin/bash

# Description: Wrapper to run R scripts

# Check the number of arguments
# if [[ -z $1 ]]
# then
#     echo "No filename for source file was po. Exiting..."
#     exit 1
# else
#     echo "$1 will be used as basename for input/output"
#     ibase=$1
# fi
# Import variables
source Env_variables/Degradome_Oliver-2022.txt 

# Activate environment
source activate ${conda_env_main}
dir=$(pwd)
# echo $dir

# Run RScripts
# echo "Running 00-Initialization.R"
# Rscript ${script_dir}/R/00-Initialization.R -d $dir -b $ibase

echo "Running 01-Annotation-shared_peak_classification.R"
Rscript ${script_dir}/R/01-Annotation-shared_peak_classification.R -d $dir

echo "Running 02-Filtering_Classification_Pooling.R"
Rscript ${script_dir}/R/02-Filtering_Classification_Pooling.R -d $dir

echo "Running 03-Drawing_Dplots.R"
Rscript ${script_dir}/R/03-Drawing_Dplots.R -d $dir

echo "Running 04-GetPeakSeq.R "
Rscript ${script_dir}/R/04-GetPeakSeq.R -d $dir

echo "Running 05-Peak_miRNA_alignment.py"
python ${script_dir}/Aux/05-Peak_miRNA_alignment.py -rd $dir

conda deactivate
echo "Running 05-Report-settings.R"
Rscript ${script_dir}/R/05-Report-settings.R -d $dir

