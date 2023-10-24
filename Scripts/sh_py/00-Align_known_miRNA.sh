#!/bin/bash

# Script to download alignments between miRNA species and their known targets
# Data comes from tarDB (Liu et. al. 2021 BMC Genomics 22 (1): 348).
# A hand-made list of species with available data is used to download,
# a file with miRNA species and their target sequence
# The above file is used to align the sequences using a global alignment and
# a custom library 'mirmap' (https://mirmap.ezlab.org/)
# Vejnar and Zdobnov, 2012 Nucleic Acids Research 40: 11673-11683

# Output: If data is available for the species under study then,
# the following files are produced.
# /Genetic_data/Fasta/<sp>_mir.fa: miRNA mature sequence
# /Genetic_data/Fasta/<sp>_tar.fa: Species target sequence
# Genetic_data/Others/Alignment_score_summary_sp_thaliana_miRNA.txt:
# Global alignment summary table
# Genetic_data/Others/DeltaG_duplex_summary_sp_thaliana_miRNA.txt:
# mirmap alignment summary table
# Genetic_data/Others/mirmap_indices_sp_thaliana_miRNA.obj
# mirmap alignments objects

# Execution: ./Scripts/sh_py/00-Align_known_miRNA.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'

# Example: /bin/bash ./Scripts/sh_py/00-Align_known_miRNA.sh Oliver-2022 Oliver-2022_vars.txt
#==================================================
# Import variables
ivars=Env_variables/Degradome_"$1".txt
if [[ ! -f "$ivars" ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh "$1" "$2"
    source "$ivars"
else
    source "$ivars"
fi

# Enable extended pattern matching
shopt -s extglob

eval "$(conda shell.bash hook)"
source activate "$conda_pydeg_map"

#Abort if get any error
set -eo pipefail

download_dir=Genetic_data/Compressed
dest_dir=Genetic_data/Fasta

# Splits species name and changes first character to uppercase
iSp=$(echo "$sp" | sed 's/_/ /' | sed -e 's/^./\U&\E/g')

# Searches for the species name and retrieves abbreviation
read -r miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($1 ~ sp) {print $2}}'\
			"$base_dir"/Env_variables/"${tarDB_list}")

#Download database from TarDB
wget http://www.biosequencing.cn/TarDB/download/${miRNAsp}.zip -O ${download_dir}/${miRNAsp}.zip
unzip -x ${download_dir}/${miRNAsp}.zip -d ${dest_dir}

#Create fasta files from top and bottom sequences
cd ${dest_dir}
awk -v sp="${miRNAsp}" '{if ($1 ~ "^" sp) {print ">",$1,$2} else {if ($1 ~ /^tar/) {gsub("U","T",$3);print $3}}}' ${miRNAsp}/${miRNAsp}.cons | sed 's/-//g' > ${miRNAsp}_tar.fa
awk -v sp="${miRNAsp}" '{if ($1 ~ "^" sp) {print ">",$1,$2} else {if ($1 ~ /^mir/) {gsub("U","T",$3);print $3}}}' ${miRNAsp}/${miRNAsp}.cons | sed 's/-//g' > ${miRNAsp}_mir.fa

# Align sequences using global alignment
python Scripts/sh_py/00-miRNA_global_alignment.py -rd $(pwd) -bn $ibase


# Align sequences using mirmap
set +u
conda deactivate
source activate "${conda_pydeg_run}"
set -u
python Scripts/sh_py/00-mirmap_miRNA_alignment.py -rd $(pwd) -bn $ibase
