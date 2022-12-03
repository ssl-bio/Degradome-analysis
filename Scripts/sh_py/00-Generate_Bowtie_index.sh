#!/bin/bash

# Generates bowtie indices for mapping to genome, transcript and to remove ribosomal RNA

# Usage: ./Scripts/sh_py/00-Generate_Bowtie_index.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'

# Example: /bin/bash ./Scripts/sh_py/00-Generate_Bowtie_index.sh Oliver-2022 Oliver-2022_vars.txt
#==================================================
#function to test/make dirs
dir_exist () {
    if [[ ! -d $1 ]]
    then
        mkdir -p $1
    fi
}

# Import variables
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh ${1} ${2}
else
    source ${ivars}
fi

#activate conda environment
eval "$(conda shell.bash hook)"
source activate ${conda_pydeg_map}

#Abort if get any error
set -eo pipefail

index_dir=Genetic_data/Index
fasta_dir=Genetic_data/Fasta
fasta_list=("${At_genome}" \
		"${At_transcript}" \
		"${At_ncRNA}")
index_path_list=(${bowtie_index_genome} ${bowtie_index_Tx} ${bowtie_index_ncRNA})

END=$((${#fasta_list[@]}-1))
for i in $(seq 0 $END)
do
    bowtie_dir=$(basename ${index_path_list[i]})
    if [ ! -d ${index_dir}/${bowtie_dir} ]
    then
	dir_exist ${index_dir}/${bowtie_dir}
	bowtie2-build ${fasta_list[i]} ${index_path_list[i]}
	cp ${At_genome} ${index_dir}/${bowtie_dir}/${bowtie_dir}.fa
	wait
    fi
done
