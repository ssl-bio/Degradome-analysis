## 02-Filtering_Classification_Pooling.R

## Description: Filters the pydeg processed output (01-pyDegradome_processed/*/*_Processed) by height of peak and classifies the resulting candidate peaks in two steps:
## Classification 1: Four categories. Based on whether the peak is found in two replicates (shared) and whether its signal to noise ratio is above certain threshold
## Classification 2: Three categories. Calculates the mean peak signal in two test replicates and compares it to the mean of the highest noise signal in control replicates.

## Output dir: output_02/02-pyDegradome_pooled
## Output files:
## -  Classification_<sample1_rep1>-<sample2_rep1>_and_<sample1_rep2>-<sample2_rep2>_<conf>_<win>_<mf>
## -  <sample1_rep1>-<sample2_rep1>_<conf>_<win>_<mf>_Annotated
## Pydegradome settings: conf: Confidence level; win: width of peak search window; mf: Multiplicative factor.

## Example: output_02/01-pyDegradome_processed/SRR10759112-SRR10759114/SRR10759112-SRR10759114_0_95_4_4_Annotated

## Libraries
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgsub))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PostPyDeg))
suppressPackageStartupMessages(library(optparse))

## Parse arguments
option_list <- list(
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Working directory", metavar="character"),
    make_option(c("-b", "--base"), type="character", default=NULL, 
                help="Project basename", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$wd)

## Import system variables
ivarF <- "Initialization_variables.RData"
ivarM <- paste0("minimal_variables_", opt$base, ".RData")
error_msg <- " was not found.\nPlease run 'Scripts/R/00-Initialization.R' providing the root directory and the project base name.\nAlternatively run '02-Degradome.sh' providing the project base name.\n"
                                        #Minimal variables
min_variables <- file.path("Env_variables", ivarM)
if (file.exists(min_variables)) {
    load(min_variables)
    ivars <- file.path(supp_data_dir,"R/Initialization_variables.RData")
    if(file.exists(ivars)) {
        load(ivars)
    } else {
        cat(ivarF, error_msg)
        break
    }
} else {
    cat(ivarM, error_msg)
    break
}
##--------------------------------------------------
                                        # Load miRNA list
mirna <- readRDS(file.path(env["supp_data_dir"],"R/miRNA"))            

                                        # Degradome BigWig
dg_bigwig_all <- list()
for (i.sample in sample_list) {
                                        #Load bigWig
    bigwig_f <- import(file.path(bigwig_dir,paste0(i.sample,".mapped_genome_f_DESeq.bw")))
    strand(bigwig_f) <- "+"

    bigwig_r <- import(file.path(bigwig_dir,paste0(i.sample,".mapped_genome_r_DESeq.bw")))
    strand(bigwig_r) <- "-"

    dg_bigwig_all[[i.sample]] <- c(bigwig_f,bigwig_r)
} 
                                        # Transcript
txdb <- loadDb(file.path(env["supp_data_dir"],"R/txdb_object"))
tr_range <- GenomicFeatures::transcripts(txdb)

                                        #Columns for index & comparison
i.cols.indx  <- c("tx_name",
                  "chr",
                  "strand",
                  "peak_start",
                  "peak_stop")

                                        # lists for summarizing peak counts
list.peak.cat1.counts <- list()
list.peak.cat2.counts <- list()
list.peak.counts <- list()
##--------------------------------------------------
for(k in seq_along(comp_pair_list)) {
    comp_list <- comp_pair_list[[k]]
    for(i.test in c("test","nc")) {
        for(i.comp in comp_list[[i.test]]) {
            for (i in seq_along(MF_list)) {
                i.MF <- MF_list[i]
                i.conf <- conf_list[i]
                i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"),c("","-"))
                i.conf_f <- gsub("\\.","_",i.conf)
                input_file <- file.path(pydeg_pooled_dir,
                                        paste("Comparison",
                                              i.comp_f,i.conf_f,"4",i.MF,sep="_"))
                out_file <- file.path(pydeg_pooled_dir,
                                      paste("Classification",i.comp_f,
                                            i.conf_f,"4",
                                            i.MF,sep="_"))
                if(!file.exists(out_file) && file.exists(input_file)) {
                    cat("Filtering and Classification 1")
                    cat("\tMF=",i.MF," conf=", i.conf, " Set=", i.test, "\n")
                    pydeg <- fread(input_file,
                                   stringsAsFactors = FALSE,
                                   header = TRUE, sep="\t")
                    ##--------------------------------------------------
                    ##Filtering parameters
                                        # Max non_peak_ratio
                    max_non_peak_ratio_thr <- 1
                    
                                        # Set Max_peak_thr   
                    max_peak_thr <- ceiling(mean(c(pydeg$max_peak_1,
                                                   pydeg$max_peak_2),
                                                 na.rm = TRUE))
                    
                                        # Peak width
                    peak_width_thr <- ceiling(mean(pydeg$width))+1
                    ##--------------------------------------------------
                    ## Filter data
                                        #Select only rows with peak:non-peak
                                        #ratio values
                    pydeg2 <- pydeg[!is.na(pydeg$max_non_peak_ratio_1) & !is.na(pydeg$max_non_peak_ratio_2),]
                    
                                        #Filter by width and max_peak count
                    ## pydeg2 <- pydeg2[pydeg2$max_peak_1 > max_peak_thr & pydeg2$max_peak_2 > max_peak_thr & pydeg2$width < peak_width_thr, ]
                    pydeg2 <- pydeg2[pydeg2$max_peak_1 > max_peak_thr & pydeg2$max_peak_2 > max_peak_thr, ]
                    
                    ##--------------------------------------------------
                    ## Classify into categories
                                        #Category 1
                    pydeg2[shared == 3 & max_non_peak_ratio_1 > 1 & max_non_peak_ratio_2 > 1,
                           category_1:= 1]

                                        #Category 2
                    pydeg2[shared == 3 & max_non_peak_ratio_1 >= 0.8 & max_non_peak_ratio_2 >= 0.8 & max_non_peak_ratio_1 <= 1 & max_non_peak_ratio_2 <= 1,
                           category_1:= 2]

                                        #Category 3
                    pydeg2[shared == 3 & max_non_peak_ratio_1 >= 0.7 & max_non_peak_ratio_2 >= 0.7 & max_non_peak_ratio_1 < 0.8 & max_non_peak_ratio_2 < 0.8,
                           category_1:= 3]

                                        #Category 4
                    pydeg2[shared %in% 1:2 & ((max_non_peak_ratio_1 >= 0.8 & max_non_peak_ratio_2 >= 1) | (max_non_peak_ratio_1 >= 1 & max_non_peak_ratio_2 >= 0.8)),
                           category_1:= 4]

                                        #Set NA in category 1 as 0
                    pydeg2[is.na(category_1),category_1:= 0 ]

                                        #------------------------------
                                        #Add miRNA info
                    sel  <- pydeg2$ID %in% mirna$ID

                    pydeg2$miRNA <- 0

                    pydeg2[sel,miRNA := 1]
                                        #------------------------------
                    cat("Classification 2")
                    i.samples <- unlist(strsplit(gsub("_and_","-",i.comp_f),
                                                 split="-"))
                    ##control samples
                    i.pairs.ctrl <- unique(i.samples[c(2,4)])
                    
                    ##test samples
                    i.pairs.test <- unique(i.samples[c(1,3)])

                    tmp <- sapply(i.samples,
                                  function(x) {
                                      as.numeric(grepl(x, names(dg_bam_list)))
                                  })
                    
                    sel_bam <- as.logical(apply(tmp,1,sum))
                    
                    bigwigs <- dg_bigwig_all[sel_bam]

                    ##Get max transcript ctrl
                    cat("Get the highest signal in the control sample\n")
                                        #Select transcripts from which to obtain
                                        #the highest read over the transcript
                    pydeg_tx <- unique(pydeg2$tx_name)
                    sel_tx <- tr_range$tx_name %in% pydeg_tx
                    gr_tr <- tr_range[sel_tx,]
                    
                    maxTxctrl_list <- maxTxctrl(tx_id=pydeg_tx,
                                                f_bigwigs=bigwigs[i.pairs.ctrl],
                                                f_pairs=i.pairs.ctrl,
                                                f_gr_tr=gr_tr,
                                                f_input_dir=maxR_dir,
                                                f_core=env$cores)
                    maxTxctrl_df <- do.call(rbind,maxTxctrl_list)
                    maxTxctrl_df <- aggregate(max_read_tx~tx_name,
                                              data=maxTxctrl_df,
                                              mean)
                    ##--------------------------------------------------
                    ##Get max peak test
                    cat("Get the peak signal in the test sample\n")
                    mydf <- dplyr::select(pydeg2,all_of(i.cols.indx))
                    maxPeakTest_list  <- maxPtest(f_df=mydf,
                                                  f_pairs=i.pairs.test,
                                                  f_bigwigs=bigwigs[i.pairs.test],
                                                  f_input_dir=maxP_dir,
                                                  f_core=env$cores)
                    maxPeakTest_df <- do.call(rbind,maxPeakTest_list)

                                        #Create an index for aggregating
                    maxPeakTest_df$indx <-  paste(maxPeakTest_df$tx_name,
                                                  maxPeakTest_df$chr,
                                                  maxPeakTest_df$strand,
                                                  maxPeakTest_df$peak_start,
                                                  maxPeakTest_df$peak_stop,sep = "_")
                    maxPeakTest_df <- aggregate(max_peak~indx,
                                                data=maxPeakTest_df,
                                                mean)
                    
                    ##-----------------------------------
                    ##Merge by TxID and calculate ratio
                    cat("Calculate ratio Peak[test]:MaxRead[ctrl]\n")
                                        #split indx to get Tx.id
                    maxPeakTest_df$tx_name <- unlist(strsplit(maxPeakTest_df$indx,"_"))[c(TRUE,FALSE,FALSE,FALSE,FALSE)]
                    pydeg.tmp <- merge(maxPeakTest_df,maxTxctrl_df,
                                       by="tx_name",
                                       all.x = TRUE)
                    pydeg.tmp$ratioPTx <- pydeg.tmp$max_peak/pydeg.tmp$max_read_tx
                                        #change column names
                    setDT(pydeg.tmp)
                    data.table::setnames(pydeg.tmp,
                                         c("tx_name",
                                           "indx",
                                           "max_peak_cpm",
                                           "max_read_tx_cpm",
                                           "ratioPTx"))
                    ##-----------------------------------
                                        #Merge with original dataframe
                    cat("Merge with input\n")
                                        #create an indx for merging
                    pydeg2 <- within(pydeg2,
                                     indx <- paste(tx_name,
                                                   chr,
                                                   strand,
                                                   peak_start,
                                                   peak_stop,
                                                   sep="_"))
                    
                    
                    pydeg2 <- merge(pydeg2,
                                    dplyr::select(pydeg.tmp,-tx_name),
                                    by="indx",
                                    all.x = TRUE)[,-1]
                    

                                        #select transcripts with more than one peak
                    dups <- data.frame(table(pydeg2$tx_name))
                    dups <- dups[dups$Freq>1,]
                    sel <- pydeg2$tx_name %in% dups$Var1
                    
                    pydeg2$MorePeaks <- as.numeric(sel)

                    ##Convert to data.table
                    pydeg2 <- as.data.table(pydeg2)
                    ##--------------------------------------------------
                    ## Classify into categories
                    
                                        #Median of ratios
                    i.RM <- median(pydeg2$ratioPTx)
                                        #Category 1
                    pydeg2[ratioPTx > 1,
                           category_2 := "A"]
                    
                                        #Category 2
                    pydeg2[ratioPTx <= 1 &
                           ratioPTx > i.RM &
                           MorePeaks>0,
                           category_2 := "B"]
                    
                                        #Set NA in category as 0
                    pydeg2[is.na(category_2),category_2 := "C" ]

                    ##--------------------------------------------------
                                        #Write filtered pydeg data frame
                    fwrite(pydeg2, out_file, quote=FALSE, sep = "\t")
                    cat("\t\tWrote",
                        paste(i.comp_f,i.conf_f,"4",i.MF,sep="_"),"\n")
                    ##--------------------------------------------------
                    ##Summarize peak numbers as they were filtered
                                        #Peak counts filtered
                    i.total.f <- nrow(pydeg2)
                    i.total <- nrow(pydeg)
                    i.settings <- paste(i.MF,i.conf,sep = "-")
                    i.cat1 <- nrow(pydeg2[pydeg2$category_1 == 1, ])
                    i.cat2 <- nrow(pydeg2[pydeg2$category_2 == "A", ])
                    i.cat1_2 <- nrow(pydeg2[pydeg2$category_1 == 1 & pydeg2$category_2 == "A", ]) 

                    l.indx <- paste0(i.comp_f,i.settings)
                    peak_counts <- data.frame(Total=i.total,
                                              Total_filtered=i.total.f,
                                              Classif_1_cat_1=i.cat1,
                                              Classif_2_cat_A=i.cat2,
                                              Cat_1_Cat_A=i.cat1_2,
                                              Comparison=i.comp_f,
                                              Settings=i.settings)
                    list.peak.counts[[l.indx]] <- peak_counts
                    
                                        #Peak counts per category 1 
                    peak_cat1_counts <- data.frame(table(pydeg2$category_1))
                    colnames(peak_cat1_counts) <- c("Category","Peaks")
                    peak_cat1_counts$Comparison <- i.comp_f
                    peak_cat1_counts$Settings <- i.settings
                    list.peak.cat1.counts[[l.indx]] <- peak_cat1_counts

                                        #Peak counts per category 2 
                    peak_cat2_counts <- data.frame(table(pydeg2$category_2))
                    colnames(peak_cat2_counts) <- c("Category","Peaks")
                    peak_cat2_counts$Comparison <- i.comp_f
                    peak_cat2_counts$Settings <- i.settings
                    list.peak.cat2.counts[[l.indx]] <- peak_cat2_counts
                }# Check for input and output
                gc()
            }# Loop over pydeg settings           
        }#i.comp
    }#Loop over test and negative control
}
##Merge peak counts and export

## Total before and after filtering
out_file <- file.path(summary_dir,"Peak_counts_BeforeAfter_Filtering")
if (!file.exists(out_file)) {
    peak_counts_i.comp <- do.call(rbind,list.peak.counts)
    peak_counts_i.comp <- peak_counts_i.comp[with(peak_counts_i.comp,order(Settings)),]


    if (!is.null(peak_counts_i.comp) && nrow(peak_counts_i.comp) > 0 ) {
        write.table(peak_counts_i.comp,
                    out_file, quote=FALSE,
                    sep="\t",row.names=FALSE)
        cat("\t\tWrote ",basename(out_file),"\n")
    }
}

## Peaks per categories
out_file <- file.path(summary_dir,"Peak_counts_classification-1")
if (!file.exists(out_file)) {
    peak_counts_cat1_i.comp <- do.call(rbind,list.peak.cat1.counts)
    peak_counts_cat1_i.comp <- peak_counts_cat1_i.comp[with(peak_counts_cat1_i.comp,order(Settings)),]



    if (!is.null(peak_counts_cat1_i.comp) && nrow(peak_counts_cat1_i.comp) > 0 ) {
        write.table(peak_counts_cat1_i.comp,
                    out_file, quote=FALSE,
                    sep="\t",row.names=FALSE)
        cat("\t\tWrote ",basename(out_file),"\n")
    }
}

## Peaks per categories
out_file <- file.path(summary_dir,"Peak_counts_classification-2")
if (!file.exists(out_file)) {
    peak_counts_cat2_i.comp <- do.call(rbind,list.peak.cat2.counts)
    peak_counts_cat2_i.comp <- peak_counts_cat2_i.comp[with(peak_counts_cat2_i.comp,order(Settings)),]



    if (!is.null(peak_counts_cat2_i.comp) && nrow(peak_counts_cat2_i.comp) > 0 ) {
        write.table(peak_counts_cat2_i.comp,
                    out_file, quote=FALSE,
                    sep="\t",row.names=FALSE)
        cat("\t\tWrote ",basename(out_file),"\n")
    }
}
##==================================================
## Pool all samples
for(i.MF in rev(MF_list)) {
    for(i.conf in rev(conf_list)) {
        if (i.MF == 4 && i.conf == 0.95 ||
            i.MF == 3 && i.conf == 0.99) {
            i.conf_f <- gsub("\\.","_",i.conf)
            pydeg_all <- data.table()
            pydeg_output_f <- file.path(pydeg_pooled_dir,
                                                    paste("Pooled",i.conf_f,"4", i.MF,sep="_"))
            for(k in seq_along(comp_pair_list)) {
                comp_list <- comp_pair_list[[k]]
                for(i.test in c("test","nc")) {
                    for(i.comp in comp_list[[i.test]]) {
                        i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"),c("","-"))
                        pydeg_input_f <-  file.path(pydeg_pooled_dir,
                                                    paste0("Classification_",
                                                           paste(i.comp_f, i.conf_f,"4", i.MF,sep="_")))
                        
                        ## if(!file.exists(pydeg_ouput_f)) {
                            if(file.exists(pydeg_input_f)) {
                                pydeg <- fread(pydeg_input_f)

                                        #Add comparison column
                                pydeg$comparison <- i.comp_f
                                
                                        #combine data simplified data
                                pydeg_all <- rbind(pydeg_all, pydeg)
                            }
                        ## }
                    }
                }
            }
                                        #Data table with less columns
            if(!is.null(pydeg_all) && nrow(pydeg_all)>0) {
                fwrite(pydeg_all,pydeg_output_f, quote = FALSE, sep='\t')
            }
        }
    }
}
