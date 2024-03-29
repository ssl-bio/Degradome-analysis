#!/bin/bash

# 01-Degradome_ch.sh
# Description: Pipeline to map reads to a genome and transcript reference, filter those to a target chromosome and then perform the identification of significant degradation fragments between pairs of samples (PyDegradome, Gaglia et al. PLOS Pathogens 2015:11)

# Execution: 01-Degradome_ch.sh <Project_name> <Variable_specification_file> The last located in 'Env_variables'
# Example: /bin/bash 01-Degradome_ch.sh Oliver-2022 Oliver-2022_vars.txt

# Output dir: output_01/
# Main output files:
# - 05-PyDegradome/t_<sample1_rep1_c_<sample2_rep1>_<conf>_<win>_<mf>

# Example: output_01/05-PyDegradome/t_SRR10759112_c_SRR10759114_0.95_4_4
#==================================================
# Check the number of arguments
if [[ -z $1 ]]
then
    echo "No basename for input/output was provided. Exiting..."
    exit 1
else
    echo "$1 will be used as basename for input/output"
    ibase=$1
fi

# Import variables
ivars=Env_variables/Degradome_"$ibase".txt
if [[ ! -f "$ivars" ]]
then
    /bin/bash Scripts/sh_py/00-Variable_setup.sh "$ibase" "$2"
    source "$ivars"
else
    source "$ivars"
fi

#Check files
ifiles=(
    "$At_genome"
    "$At_ncRNA"
    "$At_transcript"
    "$ref_gff"
    "$ref_gff_rep"
    "$ref_gtf"
    "$ref_gtf_sorted")

file_stp=0
for ifile in "${ifiles[@]}"
do
    if [ -f "$ifile" ]
    then
	file_stp=$((file_stp+1))
    fi
done
if [ "$file_stp" -lt ${#ifiles[@]} ]
then
    /bin/bash Scripts/sh_py/00-Download-genetic-data.sh "$ibase" "$ivars"
fi

# Check for bowtie indices
indices=(
    "$bowtie_index_Tx"
    "$bowtie_index_genome"
    "$bowtie_index_ncRNA")
indx_stp=0
for indx in "${indices[@]}"
do
    if [ -d "$(dirname "$indx")" ] && [ -n "$(ls -A "$(dirname "$indx")")" ]
    then
	indx_stp=$((indx_stp+1))
    fi
done
if [ "$indx_stp" -lt ${#indices[@]} ]
then
    /bin/bash Scripts/sh_py/00-Generate_Bowtie_index.sh "$ibase" "$ivars"
fi

##Step counter
stp=1

#Functions
#function to test/make dirs
dir_exist () {
    if [[ ! -d $1 ]]
    then
        mkdir -p "$1"
    fi
}

# Loop for quality control
qc_loop() { 
    for ifastq in $(ls "$1" | grep -E '\.fq$|\.fastq$|\.fastq.gz$|\.fq.gz$')
    do
	# split filename by dot
	readarray -d . -t strarr  < <(printf '%s' "$ifastq")
	iext=${strarr[-1]}

	# account for compressed files
	if [ "$iext" != "fastq" ] && [ "$iext" != "fq"  ]
	then
	    iext="${strarr[-2]}.${strarr[-1]}"
	fi
	fbase=${ifastq/.$iext/}
	outfile=$2/"$fbase"_fastqc.html
	if [ ! -f "$outfile" ]
	then
	    #no-group = Disable grouping of bases for reads >50bp.
	    fastqc --nogroup -o "$2" "$1/$ifastq"
	fi
    done
}
#==================================================
# Enable extended pattern matching
shopt -s extglob

# #activate conda environment
eval "$(conda shell.bash hook)"
source activate "$conda_pydeg_map"

#Abort if get any error
set -eo pipefail

# Change to output dir
cd "$output_dirB"

QC_dir="00-Quality_control"
Log_dir="$QC_dir/Logs"
dir_exist $Log_dir

## 01 - Fastq quality control 1
echo "$stp - Fastq quality control 1"
stp=$((stp+1))
idir_prev="01-fastq"
idir="$QC_dir/01-fastqQC_raw"
dir_exist $idir
qc_loop $idir_prev $idir

## 02 - Filter by length, trimm adapter and quality filter
echo "$stp - Trim adapter"
stp=$((stp+1))
idir_prev="01-fastq"
idir="02-fastq_Len-Trim-Q_filtered"
dir_exist "$idir"
for ifastq in $(ls "$idir_prev" | grep -E '\.fq$|\.fastq$|\.fastq.gz$|\.fq.gz$')
do
    readarray -d . -t strarr < <(printf '%s' "$ifastq")
    fbase=${strarr[0]}
    if [ ! -f "$idir/$fbase"_trimmed.fq.gz ]
    then
	trim_galore --cores "$cores" --gzip "$idir_prev/$ifastq" --output_dir "$idir"
    fi
done

## 03 - Fastq quality control 2
echo "$stp - Fastq quality control 2"
stp=$((stp+1))
idir_prev="02-fastq_Len-Trim-Q_filtered"
idir="$QC_dir/02-fastqQC_filtered"
dir_exist "$idir"
qc_loop "$idir_prev" "$idir"


## 04 - Remove rRNA
echo "$stp - Remove rRNA"
stp=$((stp+1))
idir_prev="02-fastq_Len-Trim-Q_filtered"
idir="03-fastq_No-rRNA"
dir_exist "$idir"
idir_supp="03-Supp_info"
dir_exist "$idir_supp"
for ifastq in $(ls "$idir_prev" | grep -E '\.fq$|\.fastq$|\.fastq.gz$|\.fq.gz$')
do
    seedsub=20
    outfile="$idir/${ifastq%_trimmed.fq.gz}"_No-rRNA.fastq 
    if [ ! -f "$outfile" ]
    then
	nice bowtie2  -L "$seedsub" --threads "$cores" -U "$idir_prev/$ifastq" \
             --norc \
             -N 0 --no-1mm-upfront \
             -x "$base_dir/$bowtie_index_ncRNA"  \
             --un "$outfile" \
             -S "$idir_supp/${ifastq%_trimmed.fq.gz}".sam \
             2> "$idir_supp/${ifastq%_trimmed.fq.gz}"No-rRNA_bowtie2.log
    fi
done

## 05 - Fastq quality control 3
echo "$stp - Fastq quality control 3"
stp=$((stp+1))
idir_prev="03-fastq_No-rRNA"
idir="$QC_dir/03-fastqQC_NorRNA"
dir_exist "$idir"
qc_loop "$idir_prev" "$idir"

## 06 - Map to genome
echo "$stp - Map to genome"
stp=$((stp+1))
set +u
conda deactivate
source activate "$conda_pydeg_run"
set -u

idir_prev="03-fastq_No-rRNA"
idir="04_1-bam_genomic"
dir_exist "$idir"
for ifastq in $(ls "$idir_prev" | grep -E '\.fq$|\.fastq$|\.fastq.gz$|\.fq.gz$')
do
    readarray -d . -t strarr <<< "$ifastq"
    fbase=${strarr[0]}
    outfile="$idir/${ifastq%_No-rRNA.fastq}"_mapped_genome.bam

    if [ ! -f "$outfile" ]
    then
	tophat2 --num-threads "$cores" \
	       --read-mismatches 2 \
	       --transcriptome-max-hits 1 \
	       --GTF "$base_dir/$ref_gff" \
	       --min-intron-length 25 \
	       --max-intron-length 3000  \
               --output-dir  "$idir/${ifastq%_No-rRNA.fastq}" \
	       "$base_dir/$bowtie_index_genome" \
	       "$idir_prev/$ifastq"
	
	ireads=$(cat "$idir/${ifastq%_No-rRNA.fastq}"/align_summary.txt | grep Mapped | awk '{print $3}')
	echo "$fbase", "$ireads" >> "$idir"/Log_reads.txt
	
	mv "$idir/${ifastq%_No-rRNA.fastq}"/accepted_hits.bam "$outfile"
	# samtools index "$outfile"
    fi
done

## 07 - Filter to target chromosome and convert bam to sam
echo "$stp - Filter to target chromosome and convert bam to sam"
stp=$((stp+1))

set +u
conda deactivate
source activate "$conda_pydeg_map"
set -u

idir_prev="04_1-bam_genomic"
idir_bam="04_2-bam_chromosome"
dir_exist "$idir_bam"
idir_sam="04_3-sam_chromosome"
dir_exist "$idir_sam"

for ibam in "$idir_prev"/*.bam
do
    if [ ! -f ${ibam}.bai ]
    then
	samtools index "$ibam"
    fi
    
    ifile=$(basename "$ibam")
    outfile="$idir_sam/${ifile%_mapped_genome.bam}"_mapped_chromosome.sam
    if [ ! -f "$outfile" ]
    then
	samtools view -@ "$cores" "$idir_prev/$ifile" "$ich" > "$outfile"
    fi

    outfile="$idir_bam/${ifile%_mapped_genome.bam}"_mapped_chromosome.bam
    if [ ! -f "$outfile" ]
    then
	samtools view -@ "$cores" -b "$idir_prev/$ifile" "$ich" > "$outfile"
	samtools index "$outfile"
    fi 
done

## 08 - PyDegradome
echo "$stp - PyDegradome"
stp=$((stp+1))

# Check if sorted annotation file exists
if [[ ! ${ref_gtf_sorted:+1} ]]
then
    echo "Generate auxiliary annotation files"
    # Change to input dir
    cd "$base_dir"
    /bin/bash Scripts/sh_py/00-Process_Annotation_files.sh "$ibase" "$ivars"
    source "$ivars"
    # Change to output dir
    cd "$output_dirB"
fi

set +u
conda deactivate
source activate "$conda_pydeg_run"
set -u

idir_prev="04_3-sam_chromosome"
idir="05-PyDegradome"
dir_exist "$idir"

# Calculate fasta length
read -r fasta_len < <(awk '/^>/ {seqtotal+=seqlen; seqlen=0; seq+=1; next} {seqlen += length($0)} END{print 0.01*(seqtotal+seqlen)}' "$base_dir/$At_genome")

for isettings in "${pydeg_script_settings[@]}"
do
    for icomparisons in "${pydeg_comp_list[@]}"
    do
	icomp=(${icomparisons[@]})
	iset=(${isettings[@]})

	#First comparison
	outfile="$idir"/t_${icomp[0]}_c_${icomp[1]}_${iset[0]}_4_${iset[1]}
	if [ ! -f "$outfile"_test.tab ]
	then
	    echo "First comparison with MF=${iset[0]} and SL=${iset[1]}"
	    python "$base_dir/$pydeg_script" \
		   -gtf "$base_dir/$ref_gtf_sorted" \
		   -ctrl "$idir_prev/${icomp[1]}"_mapped_chromosome.sam  \
		   -test "$idir_prev/${icomp[0]}"_mapped_chromosome.sam \
		   -iconf "${iset[0]}" \
		   -t "$fasta_len" \
		   -w 4 \
		   -imf "${iset[1]}" \
		   -o "$outfile" &
	fi
	
	#Control comparison
	outfile="$idir"/t_${icomp[1]}_c_${icomp[0]}_${iset[0]}_4_${iset[1]}
	if [ ! -f "$outfile"_test.tab ]
	then
	    echo "Oposite comparison with MF=iset[0] and SL=iset[1]"
	    python "$base_dir/$pydeg_script" \
		   -gtf "$base_dir/$ref_gtf_sorted" \
		   -ctrl "$idir_prev/${icomp[0]}"_mapped_chromosome.sam  \
		   -test "$idir_prev/${icomp[1]}"_mapped_chromosome.sam \
		   -iconf "${iset[0]}" \
		   -t "$fasta_len" \
		   -w 4 \
		   -imf "${iset[1]}" \
		   -o "$outfile"
	fi
    done
done
wait

set +u
conda deactivate
source activate "$conda_pydeg_map"
set -u

## 09 - Bam Quality check
echo "$stp - Bam Quality check"
stp=$((stp+1))

idir_prev="04_1-bam_genomic"
idir="$QC_dir/04-bam_genome_qualimap"
dir_exist "$idir"

for ibam in "$idir_prev"/*.bam
do
    ifile=$(basename "$ibam")
    out_dir="$idir/${ifile%.bam}"
    if [ ! -d "$out_dir" ]
    then
	qualimap bamqc -bam "$idir_prev/$ifile"  -gff "$base_dir/$ref_gtf_sorted" -outdir "$out_dir"
    fi
done


## 10 - Map filtered reads to transcriptome
# echo "$stp - Map filtered reads to transcript"
# stp=$((stp+1))
# idir_prev="03-fastq_No-rRNA"
# idir_log="$Log_dir/06-Log_transcript_mapping"
# dir_exist "$idir_log"
# idir="06_1-sam_transcript"
# dir_exist "$idir"
# for ifastq in $(ls "$idir_prev" | grep -E '\.fq$|\.fastq$|\.fastq.gz$|\.fq.gz$')
# do
#     outfile="$idir/${ifastq%.No-rRNA.fastq}"_mapped_transcriptome.sam
#     if [ ! -f "$outfile" ]
#     then
# 	seedsub=20
# 	nice bowtie2  -L "$seedsub" -p "$cores"  -x  "$base_dir/$bowtie_index_Tx" \
# 	     --norc \
# 	     -U "$idir_prev/$ifastq" \
# 	     -S "$outfile" 2> \
# 	     "$idir_log/$ifastq"_transcript_bowtie2.log
#     fi
# done

## 11 - Filter by target chromosome and extract uniquely mapped [Transcript]
# echo "$stp - Filter by target chromosome and extract uniquely mapped [Transcript]"
# stp=$((stp+1))
# idir_prev="06_1-sam_transcript"
# idir="06_2-sam_transcript_ch"
# dir_exist "$idir"
# for isam in "$idir_prev"/*transcriptome.sam
# do
#     ifile=$(basename "$isam")
#     outfile_tmp="$idir_prev/${ifile%_mapped_transcriptome.sam}"_sort.sam
#     outfile="$idir/$ifile"
#     if [ ! -f "$outfile" ] || [ ! -f "$outfile_tmp" ]
#     then
# 	# Filter by chromosome
# 	# Change chromosome name for Arabidopsis
# 	if [ "$ich" == "Pt" ]
# 	then
# 	    ich2=C
# 	elif [ "$ich" == "Mt" ]
# 	then
# 	    ich2=M
# 	else
# 	    ich2="$ich"
# 	fi
# 	samtools view -H "$idir_prev/$ifile" |
# 	    awk -v ch="AT${ich2}" '{if ($1 ~ /@SQ/){if($2 ~ ch){print $0}} else {print $0}}'  > "$idir_prev/${ifile%_mapped_transcriptome.sam}"_ch.sam
# 	samtools view "$idir_prev/$ifile" |
# 	    awk -v ch="AT${ich2}" '{if($3 ~ ch){print $0}}'  >> "$idir_prev/${ifile%_mapped_transcriptome.sam}"_ch.sam

# 	# Sort
# 	samtools sort -@ "$cores"  "$idir_prev/${ifile%_mapped_transcriptome.sam}"_ch.sam -O sam -o "$idir_prev/${ifile%_mapped_transcriptome.sam}"_sort.sam

# 	# Filter by XS: tag
# 	samtools view -q 10 -G 4 -G 16\
# 		 -h   "$idir_prev/${ifile%.mapped_transcriptome.sam}"_sort.sam | grep  -v "XS:" > "$outfile"
#     fi
# done

## 12 - Convert sam to bam and create bam index [Transcript]
# echo "$stp - Convert sam to bam and create bam index [Transcript]"
# stp=$((stp+1))
# idir_prev="06_2-sam_transcript_ch"
# idir="06_3-bam_transcript_ch"
# dir_exist "$idir"
# for isam in "$idir_prev"/*.sam
# do
#     ifile=$(basename "$isam")
#     outfile="$idir/${ifile%.sam}".bam
#     if [ ! -f "$outfile" ]
#     then
# 	samtools view -@ "$cores" -S -b "$idir_prev/$ifile" > \
#                  "$outfile"
# 	samtools index "$outfile"
#     fi
# done

## 13 - Bam Transcript Quality check
# echo "$stp - Bam Transcript Quality check"
# stp=$((stp+1))

# idir_prev="06_3-bam_transcript_ch"
# idir="$QC_dir/06-bam_transcript_qualimap"
# dir_exist "$idir"

# for ibam in "$idir_prev"/*.bam
# do
#     ifile=$(basename "$ibam")
#     if [ ! -f "$idir/${ifile%.bam}"/qualimapReport.html ]
#     then
# 	echo "Working on $ifile"
# 	qualimap bamqc -bam "$idir_prev/$ifile"  -outdir "$idir/${ifile%.bam}"
#     fi
# done


## 14 - Count mapped reads to features
echo "$stp - Count mapped reads to features [genome]"
stp=$((stp+1))

idir_prev="04_2-bam_chromosome"
idir="07-htseq_genomic"
dir_exist "$idir"
for ibam in "$idir_prev"/*.bam
do
    ifile=$(basename "$ibam")
    outfile="$idir/${ifile%.bam}"_representative_gene_models_HTSeq-exon.tsv
    if [ ! -f "$outfile" ]
    then
	nice htseq-count -f bam -s yes -t exon -i Parent -m intersection-strict -n "$cores" \
	     --additional-attr=exon_number "$idir_prev/$ifile" "${base_dir}/$ref_gff_rep" \
	     -c "${outfile}"
    fi
done

set +u
conda deactivate
source activate "$conda_pydeg_R"
set -u

echo "$stp - Get size factor"
stp=$((stp+1))

if [ ! -f  "$output_dirB/$idir/Size-factor.txt" ]
then
    inputfile="${base_dir}/Env_variables/minimal_variables_${ibase}.RData"
    if [ ! -f "${inputfile}" ]
    then
	Rscript "$base_dir"/Scripts/R/00-Initialization.R \
		-b "$ibase" \
		-d "$base_dir"
    fi
    
    Rscript "$base_dir"/Scripts/R/A-GetSizeFactor.R \
	    -r "$output_dir_base" \
	    -d "$output_dirB/$idir" \
	    -b "$ibase" \
	    -s "$base_dir"
fi

## 15 - Count mapped reads to features [transcript]
# echo "$stp - Count mapped reads to features [transcript]"
# set +u
# conda deactivate
# source activate "$conda_pydeg_aux"
# set -u
# stp=$((stp+1))

# idir_prev="06_3-bam_transcript_ch"
# idir="08-salmon_transcript"
# dir_exist "$idir"
# for ibam in "$idir_prev"/*.bam
# do
#     ifile=$(basename "$ibam")
#     outdir="$idir/${ifile%.bam}"
#     if [ ! -d "$outdir" ]
#     then
# 	salmon quant -p "$cores" -t "$base_dir/$At_transcript" -l U -a "$idir_prev/$ifile" -o "$outdir"
# 	mv "$outdir"/quant.sf "$idir/${ifile%.bam}".txt
#     fi
# done

# set +u
# conda deactivate
# source activate "$conda_pydeg_R"
# set -u

# if [ ! -f  "$output_dirB/$idir/Size-factor.txt" ]
# then
#     inputfile="$base_dir/Env_variables/minimal_variables_$ibase.RData"
#     if [ ! -f "$inputfile" ]
#     then
# 	Rscript "$base_dir"/Scripts/R/00-Initialization.R \
# 		-b "$ibase" \
# 		-d "$base_dir"
	
#     fi
    
#     Rscript "$base_dir"/Scripts/R/A-GetSizeFactor.R \
# 	    -r "$output_dir_base" \
# 	    -d "$output_dirB/$idir" \
# 	    -b "$ibase" \
# 	    -s "$base_dir"
# fi

## 16 - Create 5 Wig track [genome]
echo "$stp - Create 5' Wig track [genome]"
set +u
conda deactivate
source activate "$conda_pydeg_map"
set -u

stp=$((stp+1))
idir_prev="04_2-bam_chromosome"
size_factors="07-htseq_genomic/Size-factor.txt"
idir="04_4-bigwig_chromosome"
dir_exist "$idir"
for ibam in "$idir_prev"/*.bam
do
    ifile=$(basename "$ibam")

    # Get the name of the sample (no suffixes nor extensions)
    isample=$(echo $ifile | sed 's/_.*$//')
    # Raw counts
    outfile="$idir"/${ifile%_mapped_*.bam}_G_r.bw
    if [ ! -f "$outfile" ]
    then
	echo "Generating $(basename "$outfile")"
	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
                    -o "$outfile" \
                    --binSize 1 --filterRNAstrand forward  \
                    --Offset 1 --normalizeUsing None --scaleFactor '-1'
    fi

    outfile="$idir/${ifile%_mapped_*.bam}"_G_f.bw
    if [ ! -f "$outfile" ]
    then
	echo "Generating $(basename "$outfile")"
	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
                    -o "$outfile" \
                    --binSize 1 --filterRNAstrand reverse  \
                    --Offset 1 --normalizeUsing None --scaleFactor '1'
    fi

    # Scale
    # Read size factor
    read -r size_factor < <(awk -v sample="${isample}" '{if ($1 ~ sample) {print $2}}' ${size_factors})
    outfile="$idir/${ifile%_mapped_*.bam}"_G_r_DESeq.bw 
    if [ ! -f "$outfile" ]
    then
	echo "Generating $(basename "$outfile")"
	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
                    -o "$outfile"\
                    --binSize 1 --filterRNAstrand forward  \
                    --Offset 1 --normalizeUsing None --scaleFactor -"$size_factor"
    fi

    outfile="$idir/${ifile%_mapped_*.bam}"_G_f_DESeq.bw
    if [ ! -f "$outfile" ]
    then
	echo "Generating $(basename "$outfile")"
	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
                    -o "$outfile" \
                    --binSize 1 --filterRNAstrand reverse  \
                    --Offset 1 --normalizeUsing None --scaleFactor "$size_factor"
    fi
done



## 16 - Create 5' Wig track transcript...
# echo "$stp - Create 5 Wig track [transcript]"
# stp=$((stp+1))

# idir_prev="06_3-bam_transcript_ch"
# idir="06_4-bigwig_transcript"
# dir_exist "$idir"

# for ibam in "$idir_prev"/*.bam
# do
#     ifile=$(basename "$ibam")
#     isample=$(echo $ifile | sed 's/\..*$//')
#     # Raw counts
#     outfile="$idir/${ifile%.mapped_*.bam}"_Tx_r.bw
#     if [ ! -f "$outfile" ]
#     then
# 	echo "Generating $(basename "$outfile")"
# 	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
#                     -o "$outfile" \
#                     --binSize 1 --filterRNAstrand forward  \
#                     --Offset 1 --normalizeUsing None --scaleFactor '-1'
#     fi

#     outfile="$idir"/${ifile%.mapped_*.bam}_Tx_f.bw
#     if [ ! -f "$outfile" ]
#     then
# 	echo "Generating $(basename "$outfile")"
# 	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
#                     -o "$outfile" \
#                     --binSize 1 --filterRNAstrand reverse  \
#                     --Offset 1 --normalizeUsing None --scaleFactor '1'
#     fi
    
#     # Scale
#     # Read size factor
#     read -r size_factor < <(awk -v sample="${isample}" '{if ($1 ~ sample) {print $2}}' ${size_factors})
#     outfile="$idir"/${ifile%.mapped_*.bam}_Tx_r_DESeq.bw 
#     if [ ! -f "$outfile" ]
#     then
# 	echo "Generating $(basename "$outfile")"
# 	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
#                     -o "$outfile"\
#                     --binSize 1 --filterRNAstrand forward  \
#                     --Offset 1 --normalizeUsing None --scaleFactor -"$size_factor"
#     fi

#     outfile="$idir/${ifile%.mapped_*.bam}"_Tx_f_DESeq.bw
#     if [ ! -f "$outfile" ]
#     then
# 	echo "Generating $(basename "$outfile")"
# 	bamCoverage -p "$cores" -b "$idir_prev/$ifile" \
#                     -o "$outfile" \
#                     --binSize 1 --filterRNAstrand reverse  \
#                     --Offset 1 --normalizeUsing None --scaleFactor "$size_factor"
#     fi
# done
