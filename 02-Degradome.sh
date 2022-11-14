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
# Example: /bin/bash 02-Degradome.sh Oliver-2022 Oliver-2022_vars.txt

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
ivars=Env_variables/Degradome_${ibase}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh ${1} ${2}
else
    source ${ivars}
fi

# Test if miRNA variables are available
echo "Test if miRNA variables are available"
if [[ ! ${ref_miRNA:+1} ]]
then
    download_dir=Genetic_data/Compressed
    dest_dir=Genetic_data/Annotation
    miRNAsp_list=organisms.txt.gz

    # Download list of organism from mirbase
    if [ ! -f ${dest_dir}/${miRNAsp_list%.gz} ]
    then
	echo "$stp - Download list of organisms from www.mirbase.org"
	wget https://www.mirbase.org/ftp/CURRENT/${miRNAsp_list} -O ${download_dir}/${miRNAsp_list}
	gzip -dk < ${download_dir}/${miRNAsp_list} > ${dest_dir}/${miRNAsp_list%.gz}
    fi

    # Match the species under study against the list of organims
    iSp=$(echo $sp | sed 's/_/ /' | sed -e 's/^./\U&\E/g')
    read miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}'\
			 ${dest_dir}/${miRNAsp_list%.gz})

    # If the species was found download the file and add it to the list of variables
    if [[ ${miRNAsp:+1} ]]
    then
	outfile=${dest_dir}/${miRNAsp}.gff3

	# Download miRNA sequences for the target specie
	if [ ! -f ${outfile} ]
	then
	    echo "$stp - Download miRNA sequences from www.mirbase.org"
	    stp=$((stp+1))
	    iSp=$(echo $sp | sed 's/_/ /' | sed -e 's/^./\U&\E/g')
	    read miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}' ${dest_dir}/${miRNAsp_list%.gz})
	    ifile=$(basename ${outfile})
	    wget https://www.mirbase.org/ftp/CURRENT/genomes/${ifile} -O ${outfile}
	fi

	#Add the path to the downloaded file to the list of variables
	echo "ref_miRNA=$outfile" >> ${base_dir}/Env_variables/Degradome_${ibase}.txt

	#Reload list of variables
	source ${ivars}
    fi
fi

#activate conda environment
eval "$(conda shell.bash hook)"
source activate ${conda_pydeg_R}

#Abort if get any error
set -eo pipefail

dir=$(pwd)

# Run RScripts
echo "Running 00-Initialization.R"
Rscript Scripts/R/00-Initialization.R -d $dir -b $ibase


echo "Running 01-Annotation-shared_peak_classification.R"
Rscript Scripts/R/01-Annotation-shared_peak_classification.R -d $dir -b $ibase

echo "Running 02-Filtering_Classification_Pooling.R"
Rscript Scripts/R/02-Filtering_Classification_Pooling.R -d $dir -b $ibase

echo "Running 03-Drawing_Dplots.R"
Rscript Scripts/R/03-Drawing_Dplots.R -d $dir -b $ibase

echo "Running 04-GetPeakSeq.R "
Rscript Scripts/R/04-GetPeakSeq.R -d $dir -b $ibase

echo "Running 05-Peak_miRNA_alignment.py"
set +u
conda deactivate
source activate ${conda_pydeg_map}
set -u
inputfile=${supp_data_dir}/miRNA_seq/miRNA_sequences.fa
if [[ ! -f ${inputfile} ]]
then
    /bin/bash Scripts/sh_py/05-Get_miRNA_seqs.sh $dir $ibase
fi
python Scripts/sh_py/05-Peak_miRNA_alignment.py -rd $dir -bn $ibase


echo "Running 05-Report_main.R"
conda deactivate
source activate ${conda_pydeg_R}
set -u
Rscript Scripts/R/05-Report_main.R -d $dir -b $ibase
