#!/bin/bash

# Generates an annotation file sorted by start and another with only the longest isoforms 
# Usage: ./Scripts/Aux/00-Process_Annotation_files.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/Aux/00-Process_Annotation_files.sh Zhang-2021 Zhang-2021_vars.txt
#==================================================
# Import variables
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/Aux/00-Variable_setup.sh ${1} ${2}
else
    source ${ivars}
fi

#function to test/make dirs
dir_exist () {
    if [[ ! -d $1 ]]
    then
        mkdir -p $1
    fi
}

#activate conda environment
eval "$(conda shell.bash hook)"
source activate ${conda_pydeg_map}

#Abort if get any error
set -eo pipefail

# cd ${genetic_data_dir}/Annotation
# idir=$(pwd)
# echo $idir
out_file=${ref_gtf%.gtf}_sorted_awk.gtf
if [ ! -a ${out_file} ]
then
    awk '/^[^#]/ {print $0}' ${ref_gtf} | \
	awk '{if ($3 != "gene") print $0}' | \
	awk '{if ($3 != "transcript") print $0}' | \
	sort -k1,1 -k4,4n -k5,5n > ${out_file}
fi

out_file=${ref_gff%.gff3}_repr.gff3
if [ ! -a ${out_file} ]
then
    agat_sp_keep_longest_isoform.pl --gff ${ref_gff} -o ${out_file}
fi
