#!/bin/bash

# Downloads genetic data: Annotation files, Genome files and creates a bowtie index for mapping
# Execution: ./Scripts/sh_py/00-Download-genetic-data.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/sh_py/00-Download-genetic-data.sh Zhang-2021 Zhang-2021_vars.txt
#==================================================
# Import variables
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh ${1} ${2}
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

##Step counter
stp=1

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

file_base=${Sp_base}${ver}
formats=("gtf" "gff3")
for iformat in ${formats[@]}
do
    dl_file=${download_dir}/${file_base}.${iformat}.gz
    outfile=${dest_dir}/${file_base}.${iformat}
    if [ ! -f ${dl_file} ] && [ ! -f ${outfile} ]
    then
	echo "$stp - Get annotation for genome"
	stp=$((stp+1))
	wget http://ftp.ensemblgenomes.org/pub/current/plants/${iformat}/${sp}/${file_base}.${iformat}.gz -O ${dl_file}
	gzip -dk < ${download_dir}/${file_base}.${iformat}.gz > ${base_dir}/Env_variables/Degradome_${ibase}.txt
    fi
done

#miRNA
miRNAsp_list=organisms.txt.gz
if [ ! -f ${dest_dir}/${miRNAsp_list%.gz} ]
then
    echo "$stp - Download list of organisms from www.mirbase.org"
    wget https://www.mirbase.org/ftp/CURRENT/${miRNAsp_list} -O ${download_dir}/${miRNAsp_list}
    gzip -dk < ${download_dir}/${miRNAsp_list} > ${dest_dir}/${miRNAsp_list%.gz}
fi

outfile=${dest_dir}/${miRNAsp}.gff3
if [ ! -f ${outfile} ]
then
    echo "$stp - Download miRNA sequences from www.mirbase.org"
    stp=$((stp+1))
    iSp=$(echo $sp | sed 's/_/ /' | sed -e 's/^./\U&\E/g')
    read miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}' ${dest_dir}/${miRNAsp_list%.gz})
    ifile=$(basename ${outfile})
    wget https://www.mirbase.org/ftp/CURRENT/genomes/${ifile} -O ${outfile}
fi
echo "ref_miRNA==$outfile" >> Degradome_${1}.txt

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

#Fasta files
download_dir=${genetic_data_dir}/Compressed
dest_dir=${genetic_data_dir}/Fasta

# List should match the names in 
fasta_list=("${Sp_base}.dna.toplevel.fa.gz"\
       "${Sp_base}.cdna.all.fa.gz"\
       "${Sp_base}.ncrna.fa.gz")
type_list=(dna cdna ncrna)
outpath_list=("${At_genome}" "${At_transcript}" "${At_ncRNA}")

END=$((${#fasta_list[@]}-1))
for i in $(seq 0 $END)
do
    dl_file=${download_dir}/${fasta_list[i]}
    if [ ! -f ${dl_file} ] || [ ! -f ${outpath_list[i]} ]
    then
	echo "$stp - Download sequence files: genome, cDNA, ncRNA"
	stp=$((stp+1))
	wget http://ftp.ensemblgenomes.org/pub/plants/current/fasta/${sp}/${type_list[i]}/${fasta_list[i]} -O ${dl_file}
	gzip -dk < ${download_dir}/${fasta_list[i]} > ${outpath_list[i]}
    fi
done
