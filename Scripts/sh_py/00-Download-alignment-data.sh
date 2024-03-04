#!/bin/bash

# Downloads alignment data using a BioProject accession number (PRJNAxxxxxx)

# Execution: ./Scripts/sh_py/00-Download-alignment-data.sh <Project_name> <PRJNA> (Optional)
# Note: If only one argument is provided, it will look for a  variable 'prjna' with the BioProject accession number defined in a Variable file. See '00-Variable_setup.sh'

# Example 1: /bin/bash ./Scripts/sh_py/00-Download-Alignment-data.sh Oliver-2022
# Example 2: /bin/bash ./Scripts/sh_py/00-Download-Alignment-data.sh Oliver-2022 PRJNA548174
#==================================================
#function to test/make dirs
dir_exist () {
    if [[ ! -d "$1" ]]
    then
        mkdir -p "$1"
    fi
}

# Import variables
ivars=Aux_files/"$1"_vars.txt
if [[ ! -f "$ivars" ]]; then
    if [[ -z $2 ]]; then
	echo "No variable file, nor project name provided. Exiting..."
	exit 1
    else
	prjna="$2"
	lite_dir=00-lite
	fastq_dir=01-fastq
    fi
else
    source "$ivars"
    lite_dir=${output_dirB}/00-lite
    fastq_dir=${output_dirB}/01-fastq
fi

dir_exist ${lite_dir}
dir_exist ${fastq_dir}
#############################################
# Download project info with download links #
#############################################
# Download project info
project_file=Aux_files/${1}_files.csv
if [ ! -f ${project_file} ]; then
    esearch -db sra -query ${prjna} | efetch -format runinfo > ${project_file}
fi

# Get the column name for the Run and Download link
header=$(head -n 1 ${project_file})

run_col=-1
dl_col=-1

IFS=',' read -ra columns <<< "$header"
for i in "${!columns[@]}"; do
    col="${columns[i],,}"
    # Check if the column starts with 'run' or 'download'
    if [[ "$col" == run ]]; then
        run_col=$i
	run_col=$((run_col+1))
    elif [[ "$col" == download* ]]; then
        dl_col=$i
	dl_col=$((dl_col+1))
    fi
done

# Download files
while IFS= read -r line
do
    srr=$(echo $line | cut -f $run_col -d ',' | grep SRR)
    dl=$(echo $line | cut -f $dl_col -d ',')
    if [ -n "$srr" ] && [ ! -f ${out_file}  ]
    then 
        echo -e "\nDownloading $srr"
	wget ${dl} -P ${lite_dir}
    fi
done < "${project_file}"
echo -e "\nFinished downloading sequencing files"

##################################
# Extract and Compress files     #
##################################
for ilite in $(ls ${lite_dir}/*lite.1)
do
    in_file=$(basename "$ilite")
    out_file=${fastq_dir}/${in_file%*lite.1}.fastq
    if [ ! -f ${out_file}.gz ]; then
	# Extract files
        fasterq-dump -e ${cores} ${ilite} -o ${out_file}
	#Compress
        pigz -p ${cores} ${out_file} || gzip ${out_file}
    fi
done
echo -e "\nFinished extracting and compressing files"
