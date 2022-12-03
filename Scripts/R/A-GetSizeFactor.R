## A-GetSizefactor.R

## Using the read counts it calculates the size factor used for normalization. This is done using DESeq2 library. The output consist in a tab-separated file with the calculated size factor for each sample.

## Output dir: output_01/07-htseq_genomic or output_01/08-salmon_transcript
## Output files:
## - Size-factor.txt

## Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))


## Parse arguments
option_list <- list(
    make_option(c("-r", "--rd"), type="character", default=NULL, 
                help="Root directory", metavar="character"),
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Destination directory", metavar="character"),
    make_option(c("-b", "--base"), type="character", default=NULL, 
                help="Project basename", metavar="character"),
    make_option(c("-s", "--sd"), type="character", default=NULL, 
                help="Script directory", metavar="character")
) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$rd)
out_file <- file.path(opt$wd,"Size-factor.txt")
if (!file.exists(out_file)) {
    ## Import system variables
    ivarF <- "Initialization_variables.RData"
    ivarM <- paste0("minimal_variables_", opt$base, ".RData")
    error_msg <- " was not found.\nPlease run 'Scripts/R/00-Initialization.R' providing the root directory and the project base name.\nAlternatively run '02-Degradome.sh' providing the project base name.\n"

    ## Minimal set of variables
    min_variables <- file.path(opt$sd,"Env_variables", ivarM)
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

    ## Merge counts
    counts_list <- list()
    list_counts <- list.files(opt$wd, pattern ="txt|tsv")
    for (base_name in sample_list) {
        ## Get input file
        isel <- grepl(base_name,list_counts)
        ifile <- list_counts[isel]
        ## Read table
        iext <- unlist(strsplit(ifile,"\\."))
        iext <- iext[length(iext)]
        if (iext == "tsv") {
            idf <- read.table(file.path(opt$wd,ifile), header=FALSE)
            colnames(idf) <- c("ID","Count")
            idf$ID <- gsub("transcript:", "", idf$ID)
        } else {
            idf <- read.table(file.path(opt$wd,ifile), header=TRUE)
            icols <- colnames(idf)[c(1,ncol(idf))]
            idf <- dplyr::select(idf,all_of(icols))
            colnames(idf) <- c("ID","Count")
        }
        counts_list[[base_name]] <- idf
    }
    counts_data <- counts_list %>% reduce(full_join, by='ID')
    rownames(counts_data) <- counts_data[,1]
    counts_data <- counts_data[,-1]
    colnames(counts_data) <- names(counts_list)

    columnData <- idesign
    counts_data <- dplyr::select(counts_data,
                                 all_of(columnData$sample))
    DESeq.ds <- DESeqDataSetFromMatrix(countData=counts_data,
                                       colData=columnData,
                                       design=~treatment)
    DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))>0, ]
    DESeq.ds <- estimateSizeFactors(DESeq.ds)
    size.factor <- sizeFactors(DESeq.ds)

    write.table(data.frame(sample=names(size.factor), 
                           sf=1/size.factor), 
                out_file, 
                row.names=FALSE, 
                quote=FALSE, 
                col.names=FALSE)

}
