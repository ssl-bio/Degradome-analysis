
#!/bin/bash

# Generates bowtie indices for mapping to genome, transcript and to remove ribosomal RNA
# Usage: ./Scripts/Aux/00-Generate_Bowtie_index.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/Aux/00-Generate_Bowtie_index.sh Zhang-2021 Zhang-2021_vars.txt
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

index_dir=${genetic_data_dir}/Index
fasta_dir=${genetic_data_dir}/Fasta
fasta_list=("${At_genome}" \
		"${At_transcript}" \
		"${At_ncRNA}")
index_path_list=(${bowtie_index_genome} ${bowtie_index_Tx} ${bowtie_index_ncRNA})

END=$((${#fasta_list[@]}-1))
for i in $(seq 0 $END)
do
    # ifasta=${fasta_list[i]}
    bowtie_dir=$(basename ${index_path_list[i]})
    if [ ! -d ${index_dir}/${bowtie_dir} ]
    then
	dir_exist ${index_dir}/${bowtie_dir}
	bowtie2-build ${fasta_list[i]} ${index_path_list[i]}
	cp ${At_genome} ${index_dir}/${bowtie_dir}/${bowtie_dir}.fa
	wait
    fi
done
