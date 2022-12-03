#!/bin/bash

# Generates an annotation file sorted by start and another with only the longest isoforms 
# Usage: ./Scripts/Aux/00-Process_Annotation_files.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/sh_py/00-Process_Annotation_files.sh Oliver-2022 Oliver-2022_vars.txt
#==================================================
# Import variables
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/Aux/00-Variable_setup.sh ${1} ${2}
    source ${ivars}
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

#Check if annotation files exist
if [[ ! -f ${ref_gtf} ]]
then
    /bin/bash Scripts/sh_py/00-Download-genetic-data.sh ${1} ${2}
fi

outfile=${ref_gtf%.gtf}_sorted.gtf
if [ ! -a ${outfile} ]
then
    awk '/^[^#]/ {print $0}' ${ref_gtf} | \
	awk '{if ($3 != "gene") print $0}' | \
	awk '{if ($3 != "transcript") print $0}' | \
	sort -k1,1 -k4,4n -k5,5n > ${outfile}
fi
#Add the path to the sorted annotation file to the list of variables
echo "ref_gtf_sorted=$outfile" >> ${base_dir}/Env_variables/Degradome_${ibase}.txt

outfile=${ref_gff%.gff3}_repr.gff3
if [ ! -a ${outfile} ]
then
    agat_sp_keep_longest_isoform.pl --gff ${ref_gff} -o ${outfile}
fi
#Add the path to the sorted annotation file to the list of variables
echo "ref_gff_rep=$outfile" >> ${base_dir}/Env_variables/Degradome_${ibase}.txt
