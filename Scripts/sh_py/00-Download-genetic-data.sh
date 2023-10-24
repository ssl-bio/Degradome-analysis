#!/bin/bash

# Downloads genetic data: Annotation files and Genome files 

# Execution: ./Scripts/sh_py/00-Download-genetic-data.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'

# Example: /bin/bash ./Scripts/sh_py/00-Download-genetic-data.sh Oliver-2022 Oliver-2022_vars.txt
#==================================================
#function to test/make dirs
dir_exist () {
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


##Step counter
stp=1

# Create directories
directory_list=("Annotation" "Compressed" "Fasta" "Index" "Others")
for idir in "${directory_list[@]}"
do
    dir_exist Genetic_data/"$idir"
done

# Get annotation files
# Genome
download_dir=Genetic_data/Compressed
dest_dir=Genetic_data/Annotation

file_base="$Sp_base""$ver"
formats=("gtf" "gff3")
for iformat in "${formats[@]}"
do
    dl_file="$download_dir"/"$file_base"."$iformat".gz
    outfile="$dest_dir"/"$file_base"."$iformat"
    if [ ! -f "$dl_file" ] && [ ! -f "$outfile" ]
    then
	echo "$stp - Get annotation for genome"
	stp=$((stp+1))
	wget http://ftp.ensemblgenomes.org/pub/current/plants/"$iformat/$sp/$file_base"."$iformat".gz -O "$dl_file"
	gzip -dk < "$download_dir"/"$file_base"."$iformat".gz > "$outfile"
    fi
done

#Fasta files
# download_dir=Genetic_data/Compressed
dest_dir=Genetic_data/Fasta

# List should match the names in 
fasta_list=("$Sp_base.dna.toplevel.fa.gz"\
       "$Sp_base.cdna.all.fa.gz"\
       "$Sp_base.ncrna.fa.gz")
type_list=(dna cdna ncrna)
outpath_list=("${At_genome}" "${At_transcript}" "${At_ncRNA}")

END=$((${#fasta_list[@]}-1))
for i in $(seq 0 $END)
do
    dl_file="$download_dir"/"${fasta_list[i]}"
    if [ ! -f "$dl_file" ] || [ ! -f "${outpath_list[i]}" ]
    then
	echo "$stp - Download sequence files: genome, cDNA, ncRNA"
	stp=$((stp+1))
	wget http://ftp.ensemblgenomes.org/pub/plants/current/fasta/"$sp/${type_list[i]}/${fasta_list[i]}" -O "$dl_file"
	gzip -dk < "$download_dir"/"${fasta_list[i]}" > "${outpath_list[i]}"
    fi
done

# miRNA
# Splits species name and changes first character to uppercase
iSp=$(echo "$sp" | sed 's/_/ /' | sed -e 's/^./\U&\E/g')

# Searches for the species name and retrieves abbreviation
read -r miRNAsp < <(awk -v sp="$iSp" 'BEGIN {FS="\t"} ; {if ($3 ~ sp) {print $1}}'\
			"$base_dir"/Env_variables/"${mirbase_list}")

# If the species was found download the file and add it to the list of variables
if [[ ${miRNAsp:+1} ]]
then
    outfile="$dest_dir"/"$miRNAsp".gff3

    # Download miRNA sequences for the target species
    if [ ! -f "$outfile" ]
    then
	echo "$stp - Download miRNA sequences from www.mirbase.org"
	stp=$((stp+1))
	ifile=$(basename "$outfile")
	# url_root=https://www.mirbase.org/ftp/CURRENT/genomes/ # Previous
	url_root=https://www.mirbase.org/download/
	wget "${url_root}${ifile}" -O "$outfile"
	if [ $? -ne 0 ]; then
	    echo "Try to download again without checking certificate"
	    wget --no-check-certificate "${url_root}${ifile}" -O "$outfile"
	fi
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

    # Edit chromosome names annotation file and create bed file
    sed  's/^chr//g' "$ref_miRNA" | \
	sed 's/^#.*$//g' | awk 'BEGIN{FS="\t";OFS="\t"} {if ($3=="miRNA"){print $1,$4-1,$5,$9,1000,$7}}' > "$out_dir_input"/miRNA.bed

    # Get sequences from genomic fasta file
    bedtools getfasta -fi "$At_genome" -bed "$out_dir_input"/miRNA.bed -name -s -fo "$out_dir_input"/miRNA.fa

    # Format fasta header to get only the name of the miRNA sequence
    # Pipe into seqkit and remove duplicates
    awk 'BEGIN {FS=";"}{if($1 ~ /^>/) {gsub("Name=","",$3);print ">",$3} else {print $0}}' "$out_dir_input"/miRNA.fa | seqkit rmdup -s > "$out_dir_input"/miRNA_sequences.fa
fi
