#!/bin/bash

# Description: Wrapper to further analyze the classification results from PyDegradome (PyDegradome, Gaglia et al. PLOS Pathogens 2015:11) in R.
# 00-Initialization.R: Generates variables and objects used in the rest of the scripts
# 01-Annotation-shared_peak_classification.R: Annotates peaks, calculates signal:noise and identifies shared peaks
# 02-Filtering_Classification_Pooling.R: Classifies peaks in two steps
# 03-Drawing_Dplots.R: Draws decay plots of top candidates in each category
# 04-GetPeakSeq.R : Obtains genomic sequences around peaks for miRNA alignment
# 05-Peak_miRNA_alignment.py: Identifies putative miRNA targets based on global alignment. It uses a shorter sequence centered on the peak
# 05-Peak_miRNA_mirmap.py: Identifies putative miRNA targets using mirmap. It uses a longer sequence (30) and alignment can be in a region different than the peak
# 05-Report_main.R: Produces an html report summarizing the filtering and listing candidates

# Execution: 02-Degradome.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash 02-Degradome.sh Oliver-2022 Oliver-2022_vars.txt

# Check the number of arguments
if [[ -z "$1" ]]
then
    echo "No basename for input/output was provided. Exiting..."
    exit 1
else
    echo "$1 will be used as basename for input/output"
    ibase="$1"
fi

# Import variables
ivars=Env_variables/Degradome_${ibase}.txt
if [[ ! -f "${ivars}" ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh "${1}" "${2}"
else
    source "${ivars}"
fi

# Enable extended pattern matching
shopt -s extglob

# Test if miRNA variables are available
echo "Test if miRNA variables are available"
if [[ ! ${ref_miRNA:+1} ]]
then
    eval "$(conda shell.bash hook)"
    source activate "$conda_pydeg_map"
    #Abort if get any error
    set -eo pipefail
    /bin/bash Scripts/sh_py/00-Download-genetic-data.sh "$ibase" "$ivars"
fi

#activate conda environment
eval "$(conda shell.bash hook)"
source activate "${conda_pydeg_R}"

#Abort if get any error
set -eo pipefail

dir=$(pwd)

# Run RScripts
echo "Running 00-Initialization.R"
Rscript Scripts/R/00-Initialization.R -d "$dir" -b "$ibase"

echo "Running 01-Annotation-shared_peak_classification.R"
Rscript Scripts/R/01-Annotation-shared_peak_classification.R -d "$dir" -b "$ibase"

echo "Running 02-Filtering_Classification_Pooling.R"
Rscript Scripts/R/02-Filtering_Classification_Pooling.R -d "$dir" -b "$ibase"

echo "Running 03-Drawing_Dplots.R"
Rscript Scripts/R/03-Drawing_Dplots.R -d "$dir" -b "$ibase"
## Delete empty plot files (<2kb) and re-plot
plot_dir="$output_dirR"/03-Report/Dplots
pfind=$(find "$plot_dir" -name "*.pdf" -type f -size -2k)
if [ "${#pfind}" -gt 0 ]
then
    echo -e "\t Delete empty plot files (<1kb) and re-plot"
    find "$plot_dir" -name "*.pdf" -type f -size -2k -delete
    Rscript Scripts/R/03-Drawing_Dplots.R -d "$dir" -b "$ibase"
fi

echo "Running 04-GetPeakSeq.R "
Rscript Scripts/R/04-GetPeakSeq.R -d "$dir" -b "$ibase"

set +u
conda deactivate
source activate "${conda_pydeg_map}"
set -u
# Check that miRNA file exists
inputfile="${supp_data_dir}"/miRNA_seq/input/miRNA_sequences.fa
if [ ! -f ${inputfile} ] || [ ! -s ${inputfile} ]
then
    /bin/bash Scripts/sh_py/05-Get_miRNA_seqs.sh "$ibase" "$2"
fi

echo "Running 05-Peak_miRNA_alignment.py"
python Scripts/sh_py/05-Peak_miRNA_alignment.py -rd "$dir" -bn "$ibase"

echo -e "Running 05-Peak_miRNA_mirmap.py"
set +u
conda deactivate
source activate "${conda_pydeg_run}"
set -u
python Scripts/sh_py/05-Peak_miRNA_mirmap.py -rd "$dir" -bn "$ibase"

echo "Running 05-Report_main.R"
set +u
conda deactivate
source activate "${conda_pydeg_R}"
set -u
Rscript Scripts/R/05-Report_main.R -d "$dir" -b "$ibase"

set +u
conda deactivate
source activate "${conda_pydeg_map}"
set -u
echo "Running 06-Process_data.py"
python Scripts/sh_py/06-Process_data.py -rd "$dir" -bn "$ibase"
