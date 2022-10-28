#!/bin/bash

# Script to setup variables for bash and R scripts
# Usage: /bin/bash Var_setting.sh Oliver-2022 Env_variables/Degradome_vars_mint.txt 
#==================================================
# Change working dir
cd ${1}

# Import variables
source Env_variables/Degradome_Oliver-2022.txt 

# Check output dir
out_dir=${supp_data_dir}/miRNA_seq
out_dir_sub=${out_dir}/input
if [[ ! -d  ${out_dir_sub} ]]
then
mkdir -p ${out_dir_sub}
fi

# Edit chromosome names annotation file and create bed file
sed  's/^chr//g' ${ref_miRNA} | \
    sed 's/^#.*$//g' | awk 'BEGIN{FS="\t";OFS="\t"} {if ($3=="miRNA"){print $1,$4-1,$5,$9,1000,$7}}' > ${out_dir_sub}/miRNA.bed

# Get sequences from genomic fasta file
bedtools getfasta -fi ${At_genome} -bed ${out_dir_sub}/miRNA.bed -name -s -fo ${out_dir_sub}/miRNA.fa

# Format fasta header to get only the name of the miRNA sequence
# Pipe into seqkit and remove duplicates
awk 'BEGIN {FS=";"}{if($1 ~ /^>/) {gsub("Name=","",$3);print ">",$3} else {print $0}}' ${out_dir_sub}/miRNA.fa | seqkit rmdup -s > ${out_dir_sub}/miRNA_sequences.fa
