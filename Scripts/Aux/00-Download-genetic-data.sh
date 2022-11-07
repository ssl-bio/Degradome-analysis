#!/bin/bash

# Downloads genetic data: Annotation files, Genome files and creates a bowtie index for mapping
# Execution: ./Scripts/Aux/00-Download-genetic-data.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/Aux/00-Download-genetic-data.sh Zhang-2021 Zhang-2021_vars.txt
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

# Create directories
directory_list=("Annotation" "Compressed" "Fasta" "Index" "Others")
for idir in ${directory_list[@]}
do
    dir_exist ${genetic_data_dir}/${idir}
done

# Get annotation files
# Genome
download_dir=${genetic_data_dir}/Compressed
dest_dir=${genetic_data_dir}/Annotation

file_base=Arabidopsis_thaliana.TAIR10.55
formats=("gtf" "gff3")
for iformat in ${formats[@]}
do
    wget http://ftp.ensemblgenomes.org/pub/current/plants/${iformat}/arabidopsis_thaliana/${file_base}.${iformat}.gz -O ${download_dir}/${file_base}.${iformat}.gz

    gzip -dk < ${download_dir}/${file_base}.${iformat}.gz > ${dest_dir}/${file_base}.${iformat}
done

#miRNA
wget https://www.mirbase.org/ftp/CURRENT/genomes/ath.gff3 -O ${dest_dir}/ath.gff3

#Fasta files
download_dir=${genetic_data_dir}/Compressed
dest_dir=${genetic_data_dir}/Fasta

fasta_list=(Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz\
       Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz\
       Arabidopsis_thaliana.TAIR10.ncrna.fa.gz)
type_list=(dna cdna ncrna)


END=$((${#fasta_list[@]}-1))
for i in $(seq 0 $END)
do
    wget http://ftp.ensemblgenomes.org/pub/plants/current/fasta/arabidopsis_thaliana/${type_list[i]}/${fasta_list[i]} -O ${download_dir}/${fasta_list[i]}
    ifile=${fasta_list[i]%.gz}
    gzip -dk < ${download_dir}/${fasta_list[i]} > ${dest_dir}/${ifile}
done

# Others (miRNA target list)
download_dir=${genetic_data_dir}/Compressed
dest_dir=${genetic_data_dir}/Others
wget ftp://ftp.arabidopsis.org/Genes/Ath_miRNAs_Konika_Chawla_20120215.xls -O ${dest_dir}/Ath_miRNAs_Konika_Chawla_20120215.xls

#Extract sheet as csv
# ssconvert -S -O 'sheet=MIR_TARGETS' ${dest_dir}/Ath_miRNAs_Konika_Chawla_20120215.xls ${dest_dir}/%s.csv
# mv ${dest_dir}/MIR_TARGETS.csv ${dest_dir}/Ath_miRNAs_Konika_Chawla_20120215.csv
