#!/bin/bash

#Pipeline to perform a degradome analysis (Gaglia-2015, Plos pathogens 11) 
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

#Set variables
cd $(dirname $2)

VARS="`set -o posix ; set`"

source $(basename $2)

grep -vFe "$VARS" <<<"$(set -o posix ; set)" |
    grep -v ^VARS= | grep -vE "^BASH|^SHLVL" |
    sed 's/\[[0-9]\]=//g' > Degradome_${1}.txt

unset VARS


#Variables
# conda_env_main="ssl-bio"
# conda_env_py2="ssl-bio-py2"
# seq_adapters="adapters.fa"
# cores=5


#Directories
# root="/media/saul/Elements/Bioinformatics/Degradome-$ibase"
# outdir="$root/output"

conda_dir="$HOME/.local/bin/miniconda3"
conda_env_main_dir="$conda_dir/envs/$conda_env_main"
# pydeg_script=" /media/saul/Documentos/Shared/Cloned/PyDegradome"
# pydeg_script="/home/saul/Scripts/ssl-bio/Gaglia-Python-PARE-like-script"
star_extra_scripts="/media/saul/Documentos/Shared/Cloned/STAR/extras/scripts"

#Variables 2
adapter_path=$(find $conda_env_main_dir -name $seq_adapters)

#Annotation files
# ref_gff="$root/references/Annotation/Arabidopsis_thaliana.TAIR10.53.gff3"
# ref_gtf_sorted="$root/references/Annotation/Arabidopsis_thaliana.TAIR10.53_sorted_awk.gtf"
# ref_gff_representative="$root/references/Annotation/Arabidopsis_thaliana.TAIR10.53_merged_loci.gff3"

#Bowtie indexes
# rRNA_bowtie_index="$root/references/Index/Bowtie_index_At_rRNA/Bowtie_index_At_rRNA"
# bowtie_index_genome="$root/references/Index/Bowtie_index_At_dna_toplevel/Bowtie_index_At_dna_toplevel"
# bowtie_index_Tx="$root/references/Index/Bowtie_index_At_cDNA/Bowtie_index_At_cDNA"

#Star index
# star_index_genome="$root/references/Index/Star_index_At_dna_2"
# star_index_Tx="$root/references/Index/Star_index_At_cDNA"

#file's name
# ifiles=("SRR10759112" "SRR10759113" "SRR10759114" "SRR10759115")

#Suffixes for mapping
imapping=("star" "tophat")

#Suffixes for subsampling
idirsuffix=("" "_sub")

## pyDegradome settings
# pydeg_script_settings=("0.95 4" "0.99 3")

## List of comparisons
# pydeg_comp_list=("SRR10759112  SRR10759114" "SRR10759113  SRR10759115")

##Step counter
stp=1
#==================================================
#Functions
#function to test/make dirs
dir_exist () {
    if [[ ! -d $1 ]]
    then
        mkdir -p $1
    fi
}

qc_loop () {
    for ifastq in $(ls $1 | grep -v txt)
    do
	# echo $ifastq
	ifile=$(basename $ifastq)
	readarray -d . -t strarr  < <(printf '%s' "$ifile")
	iext=${strarr[-1]}

	if [ $iext != "fastq" -a  $iext != "fq"  ]
	then
	    iext="${strarr[-2]}.${strarr[-1]}"
	fi
	ibase=${ifile/.$iext/}
	outfile=$2/${ibase}_fastqc.html
	if [ ! -f $outfile ]
	then
	    #no-group = Disable grouping of bases for reads >50bp.
	    fastqc --nogroup -o $2 $1/$ifastq
	fi
    done
}
#==================================================

#activate conda environment
eval "$(conda shell.bash hook)"
source activate $conda_env_main

#Abort if get any error
set -eo pipefail

cd $output_dirB

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
echo "$stp - Filter by length, trimm adapter and quality filter"
stp=$((stp+1))
idir_prev="01-fastq"
idir="02-fastq_Len-Trim-Q_filtered"
dir_exist $idir
for ifastq in $(ls $idir_prev)
do
    readarray -d . -t strarr <<< "$ifastq"
    ibase=${strarr[0]}
    if [ ! -f ${idir}/${ibase}_trimmed.fq.gz ]
    then
	ifile=$(basename $ifastq)
	trim_galore --cores $cores $idir_prev/${ifile} --output_dir $idir
    fi
done

## 03 - Fastq quality control 2
echo "$stp - Fastq quality control 2"
stp=$((stp+1))
idir_prev="02-fastq_Len-Trim-Q_filtered"
idir="$QC_dir/02-fastqQC_filtered"
dir_exist $idir
qc_loop $idir_prev $idir


## Todo: change extension from fastq to fq and compress


## 04 - Remove rRNA
echo "$stp - Remove rRNA"
stp=$((stp+1))
idir_prev="02-fastq_Len-Trim-Q_filtered"
idir="03-fastq_No-rRNA"
dir_exist $idir
idir_supp="03-Supp_info"
dir_exist $idir_supp
for ifastq in $(ls $idir_prev | grep fq.gz)
do
    ifile=$(basename $ifastq)
    seedsub=20
    outfile=$idir/${ifile%_trimmed.fq.gz}.No-rRNA.fastq 
    if [ ! -f $outfile ]
    then
	nice bowtie2  -L $seedsub --threads $cores -U $idir_prev/$ifile \
             --norc \
             -N 0 --no-1mm-upfront \
             -x $rRNA_bowtie_index  \
             --un $outfile \
             -S $idir_supp/${ifile%_trimmed.fq.gz}.sam \
             2> $idir_supp/${ifile%_trimmed.fq.gz}No-rRNA_bowtie2.log
    fi
done

## 05 - Fastq quality control 3
echo "$stp - Fastq quality control 3"
stp=$((stp+1))
idir_prev="03-fastq_No-rRNA"
idir="$QC_dir/03-fastqQC_NorRNA"
dir_exist $idir
qc_loop $idir_prev $idir

## 06 - Map to genome [tophat]
echo "$stp - Map to genome [tophat]"
stp=$((stp+1))
set +u
conda deactivate
source activate $conda_env_py2
set -u

idir_prev="03-fastq_No-rRNA"
idir="04-sam_genomic_tophat"
dir_exist $idir
for ifastq in $(ls $idir_prev)
do
    ifile=$(basename $ifastq)
    readarray -d . -t strarr <<< "$ifile"
    ibase=${strarr[0]}
    outfile=$idir/${ifile%.No-rRNA.fastq}.mapped_genome.sam
    if [ ! -f $outfile ]
    then
	tophat --num-threads $cores --read-mismatches 2 \
	       --transcriptome-only --transcriptome-max-hits 1 \
	       --GTF $ref_gff --min-intron-length 25 \
	       --max-intron-length 3000  \
               --no-convert-bam --output-dir  $idir/${ifile%.No-rRNA.fastq} \
	       $bowtie_index_genome $idir_prev/$ifile
	ireads=$(cat $idir/${ifile%.No-rRNA.fastq}/align_summary.txt | grep Mapped | awk '{print $3}')
	echo $ibase, $ireads >> $idir/Log_reads.txt
	
	mv $idir/${ifile%.No-rRNA.fastq}/accepted_hits.sam $outfile
    fi
done


set +u
conda deactivate
source activate $conda_env_main
set -u

# 07- Map to genome [star]
echo "$stp - Map to genome [star]"
stp=$((stp+1))
idir_prev="03-fastq_No-rRNA"
idir="04-sam_genomic_star_raw"
dir_exist $idir
for ifastq in $(ls $idir_prev | grep fastq)
do
    ifile=$(basename $ifastq)
    if [ ! -f $idir/${ifile}_Aligned.out.sam ]
    then
	STAR --runThreadN $cores --genomeDir $star_index_genome \
	     --sjdbGTFfile $ref_gtf_sorted \
	     --sjdbGTFtagExonParentTranscript Parent \
	     --alignIntronMin 25 --alignIntronMax 3000 \
	     --outSAMtype SAM \
	     --outSAMstrandField intronMotif \
	     --outReadsUnmapped Fastx \
	     --readFilesIn $idir_prev/$ifile \
	     --outSAMmultNmax -1 --outMultimapperOrder Random \
	     --quantMode TranscriptomeSAM GeneCounts\
	     --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 \
	     --outFileNamePrefix $idir/${ifile}_ \
	     --sjdbOverhang 50 --seedSearchStartLmax 30
    fi
done
#Summarize read numbers
if [ ! -f $idir/Log_reads.txt ]
then
    ilogs=($(ls $idir | grep Log.final.out))
    for ilog in ${ilogs[@]}
    do
	readarray -d . -t strarr <<< "$ilog"
	ibase=${strarr[0]}
	ireads=$(cat $idir/$ilog | grep "Uniquely mapped reads number" | awk 'BEGIN {FS="|"} ; {print $2}')
	echo $ibase, $ireads >> $idir/Log_reads.txt
    done
fi

## 08 - Sort sam and add XS attribute [star]
echo "$stp - Sort sam and add XS attribute [star]"
stp=$((stp+1))
idir_prev="04-sam_genomic_star_raw"
idir="04-sam_genomic_star"
# cd $idir_prev
dir_exist $idir
#Mv summary file
if [ ! -f $idir/Log_reads.txt ]
then
    mv $idir_prev/Log_reads.txt $idir
fi

for ifile in "${isamples[@]}"
do
    outfile=${ifile}.mapped_genome.sam
    if [ ! -f $idir/$outfile ]
    then
	awk '{
          if(substr($1,1,1)=="@"){print $0}
          }' $idir_prev/${ifile}*.sam > $idir/01-header.sam

	awk '{if(NF<=16 && $NF ~ "nM") print $0}' $idir_prev/${ifile}*.sam |
	    awk -v strType=1 \
		-f $star_extra_scripts/tagXSstrandedData.awk > $idir/02-XStagadded.sam
	
	awk '{if (NF == 16 && $NF ~ "XS") {print $0}}' $idir_prev/${ifile}*.sam > \
	    $idir/03-wXStag.sam

	cd $idir
	cat 01-* 02-* 03-* > ${ifile}_XS_tmp.sam
	
	#Sort files
	samtools sort  --threads $cores --output-fmt SAM -T sample.sort -o $outfile ${ifile}_XS_tmp.sam

	#Rm temp files
	rm [0-9]*-*.sam *tmp.sam
	cd ../
    fi
done

## 09 - Subsample SAM files
echo "$stp - Subsample SAM files"
stp=$((stp+1))
for imap in "${imapping[@]}"
do
    idir_prev="04-sam_genomic_$imap"
    idir="04-sam_genomic_${imap}_sub"
    dir_exist $idir
    while read -r line;
    do
	ibase=$(echo "$line" | cut -d " " -f 1)
	iprop=$(echo "$line" | cut -d " " -f 3)
	ifile=$(ls $idir_prev | grep $ibase | grep ".sam")
	infile=${idir_prev}/${ifile}
	outfile=${idir}/${ibase}_sub.sam
	if [ -a $infile ]
	then
	    if [ ! -f $outfile ]
	    then
		picard DownsampleSam I=$infile O=$outfile P=$iprop
	    fi
	fi
    done < $idir_prev/Prop_reads.txt
done

# 09 - PyDegradome
echo "$stp - PyDegradome"
stp=$((stp+1))
set +u
conda deactivate
source activate $conda_env_py2
set -u

for iDsuffix in "${idirsuffix[@]}"
do
    for imap in "${imapping[@]}"
    do
	
	idir_prev="04-sam_genomic_${imap}${iDsuffix}"
	idir="05-pyDegradome_${imap}${iDsuffix}"
	dir_exist $idir
	for isettings in "${pydeg_script_settings[@]}"
	do
	    for icomparisons in "${pydeg_comp_list[@]}"
	    do
		icomp=(${icomparisons[@]})
		iset=(${isettings[@]})

		#First comparison
		echo "First comparison with MF=${iset[0]} and SL=${iset[1]}"
		outfile=$idir/t_${icomp[0]}_c_${icomp[1]}_${iset[0]}_4_${iset[1]}
		if [ -a $outfile ]
		then
		    echo "$outfile exists"
		else
		    python $pydeg_script $ref_gtf_sorted $idir_prev/${icomp[1]}${iDsuffix}.sam  $idir_prev/${icomp[0]}${iDsuffix}.sam \
		           ${iset[0]} 4 ${iset[1]} $outfile
		fi
		
		#Control comparison
		echo "Oposite comparison with MF=iset[0] and SL=iset[1]"
		outfile=$idir/t_${icomp[1]}_c_${icomp[0]}_${iset[0]}_4_${iset[1]}
		if [ -a $outfile ]
		then
		    echo "$outfile exists"
		else
		    python $pydeg_script $ref_gtf_sorted $idir_prev/${icomp[0]}${iDsuffix}.sam  $idir_prev/${icomp[1]}${iDsuffix}.sam \
		           ${iset[0]} 4 ${iset[1]} $outfile
		fi
	    done
	done
    done
done

set +u
conda deactivate
source activate $conda_env_main
set -u

## 11 - Convert sam to bam and create bam index
echo "$stp - Convert sam to bam and create bam index"
stp=$((stp+1))
for iDsuffix in "${idirsuffix[@]}"
do
    for imap in "${imapping[@]}"
    do
	idir_prev="04-sam_genomic_${imap}${iDsuffix}"
	idir="04-bam_genomic_$imap${iDsuffix}"
	dir_exist $idir
	for isam in $(ls $idir_prev/*.sam)
	do
	    ifile=$(basename $isam)
	    outfile=$idir/${ifile%.sam}.bam
	    if [ ! -f $outfile ]
	    then
		samtools view -@ -S -b $idir_prev/$ifile > \
			 $outfile
		samtools index $outfile
	    fi
	done
    done
done

## 11 - Bam Quality check
echo "$stp - Bam Quality check"
stp=$((stp+1))
for iDsuffix in "${idirsuffix[@]}"
do
    for imap in "${imapping[@]}"
    do
	
	idir_prev="04-bam_genomic_${imap}${iDsuffix}"
	idir="$QC_dir/04-bam_genomic_qualimap_$imap${iDsuffix}"
	dir_exist $idir

	for ibam in $(ls $idir_prev/*.bam)
	do
	    ifile=$(basename $ibam)
	    out_dir=$idir/${ifile%.bam}
	    if [ ! -d $out_dir ]
	    then
		qualimap bamqc -bam $idir_prev/$ifile  -gff $ref_gtf_sorted -outdir $out_dir
	    fi
	done
    done
done

## 12 - Move transcript alignment files (star)
# echo "$stp - Move transcript alignment files (star)"
# stp=$((stp+1))
# idir_prev="04-sam_genomic_star_raw"
# idir="06-bam_transcript_star"
# dir_exist $idir
# for ibam in $(ls $idir_prev | grep SRR.*Transcriptome)
# do
#     ifile=$(echo $ibam | sed 's/_.*$//')
#     cp $idir_prev/$ibam $idir/${ifile}.bam
#     samtools sort $idir/${ifile}.bam -o $idir/${ifile}.mapped_transcript.bam
#     samtools index $idir/${ifile}.mapped_transcriptome.bam
# done

## 13 - Map filtered reads to transcript tophat
echo "$stp - Map filtered reads to transcript tophat"
stp=$((stp+1))
idir_prev="03-fastq_No-rRNA"
idir_log="$Log_dir/06-Log_transcript_mapping"
dir_exist $idir_log
idir="06-sam_transcript_tophat"
dir_exist $idir
for ifastq in $(ls $idir_prev | grep fastq)
do
    ifile=$(basename $ifastq)
    outfile=$idir/${ifile%.No-rRNA.fastq}.mapped_transcriptome.sam
    if [ ! -f $outfile ]
    then
	seedsub=20
	nice bowtie2  -L $seedsub -p $cores  -x  $bowtie_index_Tx \
	     --norc \
	     -U $idir_prev/$ifile \
	     -S $outfile 2> \
	     $idir_log/${ifile}_transcript_bowtie2.log
    fi
done

## 14 - Extract uniquely mapped [tophat]
echo "$stp - Extract uniquely mapped [tophat]"
stp=$((stp+1))
idir_prev="06-sam_transcript_tophat"
idir="06-sam_transcript_tophat_sort"
dir_exist $idir
for isam in $(ls $idir_prev/*.sam)
do
    ifile=$(basename $isam)
    outfile=$idir/${ifile}
    if [ ! -f $outfile ]
    then
	# #Sort
	samtools sort -@ $cores  $idir_prev/$ifile -O sam -o $idir/${ifile%.mapped_transcriptome.sam}_sort.sam

	# #Filter by XS: tag
	samtools view -q 10 -G 4 -G 16\
		 -h   $idir/${ifile%.sam}_sort.sam | grep  -v "XS:" > $outfile

    fi
done

## 15 - Convert sam to bam and create bam index [Tophat]
echo "$stp - Convert sam to bam and create bam index [Tophat]"
stp=$((stp+1))
idir_prev="06-sam_transcript_tophat_sort"
idir="06-bam_transcript_tophat"
dir_exist $idir
for isam in $(ls $idir_prev | grep "mapped_transcriptome.sam")
do
    ifile=$(basename $isam)
    outfile=$idir/${ifile%.sam}.bam
    if [ ! -f $outfile ]
    then
	samtools view -@ 8 -S -b $idir_prev/$ifile > \
                 $outfile
	samtools index $outfile
    fi
done

## 16 - Map to transcriptome [star]
echo "$stp - Map to transcriptome [star]"
stp=$((stp+1))
idir_prev="03-fastq_No-rRNA"
idir="06-sam_transcript_star_raw"
dir_exist $idir
for ifastq in $(ls $idir_prev | grep fastq)
do
    ifile=$(basename $ifastq)
    if [ ! -f $idir/${ifile}_Aligned.out.sam ]
    then
	STAR --runThreadN $cores --genomeDir $star_index_Tx \
	     --sjdbGTFtagExonParentTranscript Parent \
	     --alignIntronMin 25 --alignIntronMax 3000 \
	     --outSAMtype SAM \
	     --outSAMstrandField intronMotif \
	     --outReadsUnmapped Fastx \
	     --readFilesIn $idir_prev/$ifile \
	     --outSAMmultNmax -1 --outMultimapperOrder Random \
	     --outFilterMultimapNmax 1 --outFilterMismatchNmax 0 \
	     --outFileNamePrefix $idir/${ifile}_ \
	     --sjdbOverhang 50 --seedSearchStartLmax 30
    fi
done

## 17 - Convert sam to bam and create bam index [star]
echo "$stp - Convert sam to bam and create bam index [star]"
stp=$((stp+1))
idir_prev="06-sam_transcript_star_raw"
idir="06-bam_transcript_star"
dir_exist $idir
for isam in $(ls $idir_prev | grep ".sam")
do
    ifile=$(basename $isam)
    readarray -d . -t strarr <<< "$ifile"
    ibase=${strarr[0]}
    outfile=$idir/${ibase}.mapped_transcriptome.bam
    if [ ! -f $outfile ]
    then
	echo "Generating $outfile"
	# Sort sam files and move to another dir
	samtools sort -@ $cores  $idir_prev/$ifile -O sam \
		 -o $idir/${ifile%.sam}_sort.sam
	
	# Extract uniquely mapped
	samtools view -q 10 -G 4 -G 16\
		 -h   $idir/${ifile%.sam}_sort.sam | \
	    grep  -v "XS:" > ${outfile%.bam}.sam

	#Sort files
	samtools sort  --threads $cores --output-fmt SAM -T sample.sort -o ${outfile%.bam}_sort.sam ${outfile%.bam}.sam
	
	# Convert to bam
	samtools view -@ 8 -S -b ${outfile%.bam}_sort.sam > \
		 $outfile

	# Create index
	samtools index $outfile
    fi
done

## 18 - Bam Transcript Quality check #Fails for star
echo "$stp - Bam Transcript Quality check"
stp=$((stp+1))
for imap in "${imapping[@]}"
do
    idir_prev="06-bam_transcript_$imap"
    idir="$QC_dir/06-bam_transcript_qualimap_$imap"
    dir_exist $idir

    for ibam in $(ls $idir_prev | grep -E ".mapped_transcriptome.bam$")
    do
	if [ ! -f $idir/${ibam%.bam}/qualimapReport.html ]
	then
	    # ifile=$(basename $ibam)
	    echo "Working on $ifile"
	    qualimap bamqc -bam $idir_prev/$ibam  -outdir $idir/${ibam%.bam}
	fi
    done
done

## 19 - Count mapped reads to features
echo "$stp - Count mapped reads to features"
stp=$((stp+1))

for iDsuffix in "${idirsuffix[@]}"
do
    for imap in "${imapping[@]}"
    do
	idir_prev="04-bam_genomic_$imap${iDsuffix}"
	idir="07-htseq_genomic_$imap${iDsuffix}"
	# isuffix=
	dir_exist $idir
	for ibam in $(ls $idir_prev/*.bam)
	do
	    ifile=$(basename $ibam)
	    outfile=$idir/${ifile%.bam}_representative_gene_models_HTSeq-exon.txt
	    if [ ! -f $outfile ]
	    then
		nice htseq-count -f bam -s yes -t exon -i Parent -m intersection-strict -n $cores --additional-attr=exon_number $idir_prev/$ifile $ref_gff_representative > $idir/${ifile%.bam}_representative_gene_models_HTSeq-exon.txt &
	    fi
	done
	Rscript ${script_dir}/R/A-GetSizeFactor.R -r ${output_dir_base} \
		-d ${output_dirB}/${idir_prev}
	
    done
done

## 20 - Create 5 Wig track
echo "$stp - Create 5' Wig track"
stp=$((stp+1))
for iDsuffix in "${idirsuffix[@]}"
do
    for imap in "${imapping[@]}"
    do

	idir_prev="04-bam_genomic_${imap}${iDsuffix}"
	size_factors="07-htseq_genomic_$imap${iDsuffix}/Size-factor.txt"
	idir="04-bigwig_genomic_${imap}${iDsuffix}"
	dir_exist $idir
	for ibam in $(ls $idir_prev/*.bam)
	do
	    ifile=$(basename $ibam)
	    isample=$(echo $ifile | sed 's/\..*$//g')
	    # Raw counts
	    outfile=$idir/${ifile%.bam}_r.bw
	    if [ ! -f $outfile ]
	    then
		bamCoverage -p $cores -b $idir_prev/$ifile \
                            -o ${outfile} \
                            --binSize 1 --filterRNAstrand forward  \
                            --Offset 1 --normalizeUsing None --scaleFactor '-1'
	    fi

	    outfile=$idir/${ifile%.bam}_f.bw
	    if [ ! -f $outfile ]
	    then
		bamCoverage -p $cores -b $idir_prev/$ifile \
                            -o $outfile \
                            --binSize 1 --filterRNAstrand reverse  \
                            --Offset 1 --normalizeUsing None --scaleFactor '1'
	    fi
	    
	    # Scale
	    # Read size factor
	    read size_factor < <(awk -v sample="${isample}" '{if ($1 ~ sample) {print $2}}' ${size_factors})
	    outfile=$idir/${ifile%.bam}_r_DESeq.bw 
	    if [ ! -f $outfile ]
	    then
		bamCoverage -p $cores -b $idir_prev/$ifile \
                            -o $outfile\
                            --binSize 1 --filterRNAstrand forward  \
                            --Offset 1 --normalizeUsing None --scaleFactor -${size_factor}
	    fi

	    outfile=${ifile%.bam}_f_DESeq.bw
	    if [ ! -f $outfile ]
	    then
		bamCoverage -p $cores -b $idir_prev/$ifile \
                            -o $idir/$outfile \
                            --binSize 1 --filterRNAstrand reverse  \
                            --Offset 1 --normalizeUsing None --scaleFactor ${size_factor}
	    fi
	done
    done
done


## 21 - Create 5' Wig track transcript...
echo "$stp - Create 5 Wig track transcript"
stp=$((stp+1))
for imap in "${imapping[@]}"
do
    
    idir_prev="06-bam_transcript_$imap"
    idir="06-bigwig_transcript_$imap"
    dir_exist $idir

    for ibam in $(ls $idir_prev/*)
    do
	ifile=$(basename $ibam)
	## Raw counts
	echo "\t- Raw counts"
	if [ ! -f $idir/${ifile}_uni_r.bw ]
	then
	    bamCoverage --numberOfProcessors $cores --bam $idir_prev/$ifile \
			--outFileName $idir/${ifile}_uni_r.bw \
			--binSize 1 --filterRNAstrand forward  \
			--Offset 1 --normalizeUsing None --scaleFactor '-1'
	fi

	if [ ! -f $idir/${ifile}_uni_f.bw ]
	then
	    bamCoverage --numberOfProcessors $cores --bam $idir_prev/$ifile \
			--outFileName $idir/${ifile}_uni_f.bw \
			--binSize 1 --filterRNAstrand reverse  \
			--Offset 1 --normalizeUsing None --scaleFactor '1'
	fi

	## Normalized CPM
	echo "\t- Normalized CPM"
	if [ ! -f $idir/${ifile}_uni_r_CPM.bw ]
	then
	    bamCoverage --numberOfProcessors $cores --bam $idir_prev/$ifile \
			--outFileName $idir/${ifile}_uni_r_CPM.bw \
			--binSize 1 --filterRNAstrand forward  \
			--Offset 1 --normalizeUsing CPM --scaleFactor '-1'
	fi
	
	if [ ! -f $idir/${ifile}_uni_f_CPM.bw ]
	then
	    bamCoverage --numberOfProcessors $cores --bam $idir_prev/$ifile \
			--outFileName $idir/${ifile}_uni_f_CPM.bw \
			--binSize 1 --filterRNAstrand reverse  \
			--Offset 1 --normalizeUsing CPM --scaleFactor '1'
	fi
    done
done
