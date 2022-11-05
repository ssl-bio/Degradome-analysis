#!/bin/bash


#Pipeline to perform a degradome analysis
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
ivars=Env_variables/Degradome_${1}.txt
if [[ ! -f ${ivars} ]]
then
    /bin/bash Scripts/Aux/00-Variable_setup.sh ${1} Degradome_vars_mint.txt
else
    source ${ivars}
fi



# conda_dir="$HOME/.local/bin/miniconda3"
# conda_env_main_dir="$conda_dir/envs/$conda_env_main"

##Step counter
stp=1

#Functions
#function to test/make dirs
dir_exist () {
    if [[ ! -d $1 ]]
    then
        mkdir -p $1
    fi
}

# Loop for quality control
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

# Change to output dir
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
echo "$stp - Trim adapter"
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


# 07 - PyDegradome
echo "$stp - PyDegradome"
stp=$((stp+1))
set +u
conda deactivate
source activate $conda_env_py2
set -u

idir_prev="04-sam_genomic_tophat"
idir="05-pyDegradome_tophat"
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


set +u
conda deactivate
source activate $conda_env_main
set -u

## 08 - Convert sam to bam and create bam index
echo "$stp - Convert sam to bam and create bam index"
stp=$((stp+1))

idir_prev="04-sam_genomic_tophat"
idir="04-bam_genomic_tophat"
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


## 09 - Bam Quality check
echo "$stp - Bam Quality check"
stp=$((stp+1))

idir_prev="04-bam_genomic_tophat"
idir="$QC_dir/04-bam_genomic_qualimap_tophat"
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


## 10 - Map filtered reads to transcriptome tophat
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

## 11 - Extract uniquely mapped [tophat]
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

## 12 - Convert sam to bam and create bam index [Tophat]
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


## 13 - Bam Transcript Quality check
echo "$stp - Bam Transcript Quality check"
stp=$((stp+1))

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


## 14 - Count mapped reads to features
echo "$stp - Count mapped reads to features [genome]"
stp=$((stp+1))

idir_prev="04-bam_genomic_tophat"
idir="07-htseq_genomic_tophat"
# isuffix=
dir_exist $idir
for ibam in $(ls $idir_prev/*.bam)
do
    ifile=$(basename $ibam)
    outfile=$idir/${ifile%.bam}_representative_gene_models_HTSeq-exon.txt
    if [ ! -f $outfile ]
    then
	nice htseq-count -f bam -s yes -t exon -i Parent -m intersection-strict -n $cores \
	     --additional-attr=exon_number $idir_prev/$ifile $ref_gff_representative \
	     -c ${outfile}
    fi
done
Rscript ${script_dir}/R/A-GetSizeFactor.R \
	# -r ${output_dir_base} \
	-r ${output_dir_base} \
	-d ${output_dirB}/${idir_prev}
	
## 15 - Count mapped reads to features
echo "$stp - Count mapped reads to features [transcript]"
stp=$((stp+1))

idir_prev="06-bam_transcript_tophat"
idir="08-salmon_transcript_tophat"
dir_exist $idir
for ibam in $(ls $idir_prev/*.bam)
do
    ifile=$(basename $ibam)
    outdir=$idir/${ifile%.bam}
    if [ ! -f $outfile ]
    then
	salmon quant -p ${cores} -t ${At_transcript} -l U -a $idir_prev/$ifile -o ${outfile}
	mv $outdir/quant.sf $idir/${ifile%.bam}.txt
    fi
done
Rscript ${script_dir}/R/A-GetSizeFactor.R \
	# -r ${output_dir_base} \
	-r ${output_dir_base} \
	-d ${output_dirB}/${idir}


## 16 - Create 5 Wig track
echo "$stp - Create 5' Wig track"
stp=$((stp+1))
idir_prev="04-bam_genomic_tophat"
size_factors="07-htseq_genomic_tophat/Size-factor.txt"
idir="04-bigwig_genomic_tophat"
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



## 16 - Create 5' Wig track transcript...
echo "$stp - Create 5 Wig track transcript"
stp=$((stp+1))


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
