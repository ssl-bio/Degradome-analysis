## A-GetSizefactor.R

## Using the read counts it calculates the size factor used for normalization

## Libraries
library(DESeq2)
library(tidyverse)
library(optparse)


## Parse arguments
option_list <- list(
    make_option(c("-r", "--rd"), type="character", default=NULL, 
                help="Root directory", metavar="character"),
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Destination directory", metavar="character")
    ## make_option(c("-s", "--sd"), type="character", default=NULL, 
    ##             help="Destination directory", metavar="character")
) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$rd)
out_file <- file.path(opt$wd,"Size-factor.txt")
if (!file.exists(out_file)) {
                                        #Minimal variables
    load("Env_variables/mininal_variables.RData")

    ivars <- file.path(supp_data_dir,"R/Initialization_variables.RData")
    if(file.exists(ivars)) {
        load(ivars)
    } else {
        source("Scripts/R/00-Initialization.R")
    }

    ## Merge counts
    counts_list <- list()
    for (base_name in sample_list) {
        ## Get input file
        isel <- grepl(base_name,list.files(opt$wd))
        ifile <- list.files(opt$wd)[isel]
        cat(file.path(opt$wd,ifile),"\n")
        ## Read table
        idf <- read.table(file.path(opt$wd,ifile), header=FALSE)
        colnames(idf) <- c("ID","Count")
        counts_list[[base_name]] <- idf
    }
    counts_data <- counts_list %>% reduce(full_join, by='ID')
    rownames(counts_data) <- gsub("transcript:", "", counts_data[,1])
    counts_data <- counts_data[,-1]
    colnames(counts_data) <- names(counts_list)

    columnData <- idesign
    counts_data <- dplyr::select(counts_data,
                                 all_of(columnData$sample))
    DESeq.ds <- DESeqDataSetFromMatrix(countData=counts_data,
                                       colData=columnData,
                                       design=~treatment)
    DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))>0, ]
                                        #colSums(counts(DESeq.ds))
    DESeq.ds <- estimateSizeFactors(DESeq.ds)
    size.factor <- sizeFactors(DESeq.ds)

    write.table(data.frame(sample=names(size.factor), 
                           sf=1/size.factor), 
                out_file, 
                row.names=FALSE, 
                quote=FALSE, 
                col.names=FALSE)

}