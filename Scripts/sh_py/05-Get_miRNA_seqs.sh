#!/bin/bash

# Script to search and download miRNA annotation file from from www.mirbase.org
# then, using the annotation file it produces a fasta file with miRNA sequences.
# These sequences will be used for aligning peak sequences in order to find potential
# targets

# Usage: /bin/bash Scripts/sh_py/05-Get_miRNA_seqs.sh Oliver-2022 Oliver-2022_vars.txt 
#==================================================
#function to test/make dirs
dir_exist() {
    if [[ ! -d "$1" ]]
    then
        mkdir -p "$1"
    fi
}

# Import variables
ivars=Env_variables/Degradome_"$1".txt
if [[ ! -f "$ivars" ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh "$1" "$2"
    source "$ivars"
else
    source "$ivars"
fi

# miRNA
# Download list of organisms from www.mirbase.org

miRNAsp_list=organisms.txt.gz
download_dir=Genetic_data/Compressed
dest_dir=Genetic_data/Annotation

if [ ! -f "$dest_dir/$miRNAsp_list%.gz" ]
then
    echo "$stp - Download list of organisms from www.mirbase.org"
    wget https://www.mirbase.org/ftp/CURRENT/"$miRNAsp_list" -O "$download_dir/$miRNAsp_list"
    gzip -dk < "$download_dir/$miRNAsp_list" > "$dest_dir/${miRNAsp_list%.gz}"
fi

iSp=$(echo "$sp" | sed 's/_/ /' | sed -e 's/^./\U&\E/g')
read -r miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}'\
		     "$dest_dir/${miRNAsp_list%.gz}")

# If the species was found, download the file and add it to the list of variables
if [[ ${miRNAsp:+1} ]]
then
    outfile="$dest_dir/$miRNAsp".gff3

    # Download miRNA sequences for the target species
    if [ ! -f "$outfile" ]
    then
	echo "$stp - Download miRNA sequences from www.mirbase.org"
	stp=$((stp+1))
	ifile=$(basename "$outfile")
	wget https://www.mirbase.org/ftp/CURRENT/genomes/"$ifile" -O "$outfile"
    fi

    #Add the path to the downloaded file to the list of variables
    echo "ref_miRNA=$outfile" >> "$base_dir"/Env_variables/Degradome_"$ibase".txt

    #Reload list of variables
    source "$ivars"

    # Check output dir
    out_dir="$supp_data_dir"/miRNA_seq
    out_dir_input="$out_dir"/input
    out_dir_output="$out_dir"/output

    dir_exist "$out_dir_input"
    dir_exist "$out_dir_output"

    # Edit chromosome names annotation file and create bed file for 'miRNA' seqs
    sed  's/^chr//g' "$ref_miRNA" | \
	sed 's/^#.*$//g' | awk 'BEGIN{FS="\t";OFS="\t"} {if ($3=="miRNA"){print $1,$4-1,$5,$9,1000,$7}}' > "$out_dir_input"/miRNA.bed

    # Get sequences from genomic fasta file
    bedtools getfasta -fi "$At_genome" -bed "$out_dir_input"/miRNA.bed -name -s -fo "$out_dir_input"/miRNA.fa

    # Format fasta header to get only the name of the miRNA sequence
    # Pipe into seqkit and remove duplicates
    awk 'BEGIN {FS=";"}{if($1 ~ /^>/) {gsub("Name=","",$3);print ">",$3} else {print $0}}' "$out_dir_input"/miRNA.fa | seqkit rmdup -s > "$out_dir_input"/miRNA_sequences.fa
fi

