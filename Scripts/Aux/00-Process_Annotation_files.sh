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

if [ ! -a ${ref_gtf_sorted} ]
then
    echo "Sorting annotation (gtf) file"
    awk '/^[^#]/ {print $0}' ${ref_gtf} | \
	awk '{if ($3 != "gene") print $0}' | \
	awk '{if ($3 != "transcript") print $0}' | \
	sort -k1,1 -k4,4n -k5,5n > ${ref_gtf_sorted}
fi

if [ ! -a $ref_gff_rep} ]
then
    echo "Generating annotation file (gff3) of representative (longest) isoforms"
    agat_sp_keep_longest_isoform.pl --gff ${ref_gff} -o ${ref_gff_rep}
    mv ${base_dir}/*.agat.log ${base_dir}/Genetic_data/Annotation/
fi
