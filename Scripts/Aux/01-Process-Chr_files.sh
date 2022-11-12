#!/bin/bash

# Downloads and generates annotation and sequence files to analyze data on a single chromosome
# Usage: ./Scripts/Aux/01-Process-Chr_files.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash ./Scripts/Aux/01-Process-Chr_files.sh Zhang-2021_4 Zhang-2021_vars_ch.txt

#Output files:
# -  
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

##Step counter
stp=1

# Create directories
directory_list=("Annotation" "Compressed" "Fasta" "Index" "Others")
for idir in ${directory_list[@]}
do
    dir_exist ${genetic_data_dir}/${idir}
done

# Get annotation files
# Chromosomes
download_dir=${genetic_data_dir}/Compressed
dest_dir=${genetic_data_dir}/Annotation

file_base=${Sp_base}${ver}
format="gff3"
outfile=${ref_gff}
dl_file=${download_dir}/${file_base}${i_chr}.${format}.gz
if [ ! -f ${dl_file} ] || [ ! -f ${outfile} ]
then
    echo "$stp - Get annotation for chromosome"
    stp=$((stp+1))
    wget http://ftp.ensemblgenomes.org/pub/current/plants/${format}/${sp}/${file_base}${chr_pref}${i_chr}.${format}.gz -O ${dl_file}
    gzip -dk < ${dl_file} > ${outfile}
fi

# Convert GFF3 to GTF
if [ ! -f ${ref_gtf} ]
then
    echo "$stp - Convert GFF3 to GTF"
    stp=$((stp+1))
    agat_convert_sp_gff2gtf.pl --gff ${ref_gff} -o ${ref_gtf}
    mv ${base_dir}/*.agat.log ${base_dir}/Genetic_data/Annotation/
fi

# Sort annotation file
if [ ! -f ${ref_gtf_sorted} ]
then
    echo "$stp - Sort annotation file"
    stp=$((stp+1))
    awk '/^[^#]/ {print $0}' ${ref_gtf} | \
	awk '{if ($3 != "gene") print $0}' | \
	awk '{if ($3 != "transcript") print $0}' | \
	sort -k1,1 -k4,4n -k5,5n > ${ref_gtf_sorted}
fi

# Keep longest isoform (representative)
if [ ! -f ${ref_gff_rep} ]
then
    echo "$stp - Keep longest isoform"
    stp=$((stp+1))
    agat_sp_keep_longest_isoform.pl --gff ${ref_gff} -o ${ref_gff_rep}
    mv ${base_dir}/*.agat.log ${base_dir}/Genetic_data/Annotation/
fi

#miRNA
miRNAsp_list=organisms.txt.gz
if [ ! -f ${dest_dir}/${miRNAsp_list%.gz} ]
then
    echo "$stp - Download list of organisms from www.mirbase.org"
    wget https://www.mirbase.org/ftp/CURRENT/${miRNAsp_list} -O ${download_dir}/${miRNAsp_list}
    gzip -dk < ${download_dir}/${miRNAsp_list} > ${dest_dir}/${miRNAsp_list%.gz}
    
fi


echo "$stp - Download miRNA sequences from www.mirbase.org"
stp=$((stp+1))
iSp=$(echo $sp | sed 's/_/ /' | sed -e 's/^./\U&\E/g')
read miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}' ${dest_dir}/${miRNAsp_list%.gz})
outfile=${dest_dir}/${miRNAsp}.gff3
if [ ! -f ${outfile} ]
then
    ifile=$(basename ${outfile})
    wget https://www.mirbase.org/ftp/CURRENT/genomes/${ifile} -O ${outfile}
fi
echo "ref_miRNA==$outfile" >> Degradome_${1}.txt


# Download sequence files: genome, cDNA, ncRNA
fasta_list=("${Sp_base}.dna.toplevel.fa.gz"\
		"${Sp_base}.cdna.all.fa.gz"\
		"${Sp_base}.ncrna.fa.gz")
type_list=(dna cdna ncrna)
outpath_list=("${At_genome_aux}" "${At_transcript_aux}" "${At_ncRNA}")

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

# Subset genome for chromosome
if [ ! -f ${At_genome} ]
then
    echo "$stp - Subset genome for target chromosome"
    stp=$((stp+1))
       bioawk -v ichr=":$ich:" -c fastx '{if ($0 ~ ichr) {print ">"$1,$4; print $2}}' ${At_genome_aux} > ${At_genome}
fi

# Subset cDNA for chromosome
if [ ! -f ${At_transcript} ]
then
    echo "$stp - Subset cDNA for target chromosome"
    stp=$((stp+1))
       bioawk -v ichr=":$ich:" -c fastx '{if ($0 ~ ichr) {print ">"$1,$4; print $2}}' ${At_transcript_aux} > ${At_transcript}
fi


