## Libraries
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(optparse))

## Parse arguments
option_list <- list(
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Working directory", metavar="character")
) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$wd)

## Import system variables
                                        #Minimal variables
load("Env_variables/mininal_variables.RData")

ivars <- file.path(supp_data_dir,"R/Initialization_variables.RData")
if(file.exists(ivars)) {
    load(ivars)
} else {
    source("Scripts/R/00-Initialization.R")
}
##-------------------------------------------------------------------------

EDB <- EnsDb(DB)

i.fasta <- readDNAStringSet(env[["At_transcript"]])
names(i.fasta) <-  vapply(strsplit(names(i.fasta)," "), `[`, 1, FUN.VALUE=character(1))
fasta_names <- names(i.fasta)
peak_region <- 30
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    i.conf_f <- gsub("\\.","_",i.conf)

                                        #Check if file exists
    pydeg_input_f <- file.path(pydeg_pooled_dir,paste("Pooled",
                                                      i.conf_f, "4", i.MF, sep="_"))
    if (file.exists(pydeg_input_f)) {
        ## df.f <- paste(i.comp,i.conf,"4",i.MF,"Motifs_category_1", sep = "_")
        fasta.f <- paste0(paste("PeakRegioncDNA_category_1",i.conf_f,"4",i.MF,
                                sep = "_"),".fa")
        
                                        #Read pydeg input file
        pydeg_all <- fread(pydeg_input_f)
        sel <- pydeg_all$category_1=="1" & pydeg_all$category_2=="A"
        pydeg_sub <- pydeg_all[sel,]
        pydeg_sub$indx <- paste0(pydeg_sub$tx_name,
                                 pydeg_sub$peak_start,
                                 pydeg_sub$peak_stop)
        
        gnm <- GRanges(pydeg_sub$chr,
                       IRanges(start=pydeg_sub$peak_start,
                               end=pydeg_sub$peak_stop),
                       strand=pydeg_sub$strand)
        gnm_tx <- genomeToTranscript(gnm, EDB)  

        tx_names <- pydeg_sub$tx_name

        seq.df <- data.frame()
        for (i in seq_along(gnm_tx)) {
            itx <- tx_names[i]
            idf <- data.frame(gnm_tx[[i]])
            seltx <- idf$tx_id %in% itx
            selfa <- fasta_names %in% itx
            if (sum(seltx)==1 && sum(selfa)==1) {
                                        # Get sequence
                ifasta <- i.fasta[selfa]
                
                                        # Get coordinates around peak
                idf_sub <- idf[seltx,]
                iwdth <- idf[seltx,"width"]
                if(iwdth%%2 == 0) {
                    i.off <- (peak_region-iwdth)/2
                } else {
                    i.off <- (peak_region-iwdth-1)/2
                }
                i.start <- idf_sub$start - i.off
                i.end <-  idf_sub$end + i.off
                if(i.start < 1) {
                    i.start <- 1
                }
                if(i.end > nchar(ifasta)) {
                    i.end <- nchar(ifasta)
                }
                                        # Get sequence around peak
                i.seq <- subseq(ifasta, start=i.start, end=i.end)
            } else {
                idf_sub <- idf[1,]
                i.seq <- NA
                i.start <-NA
                i.end <- NA
            }

            idf_sub$TxP_start <- i.start
            idf_sub$TxP_end <- i.end
            idf_sub$Pseq <- as.character(i.seq)
            idf_sub$indx <- paste0(idf_sub$tx_id,idf_sub$seq_start,idf_sub$seq_end)
            ## idf <- data.frame(tx_name=itx,
            ##                   TxP_start=i.start,
            ##                   TxP_end=i.end,
            ##                   Pseq=i.seq)
            seq.df <- rbind(seq.df,idf_sub)
        }
        seq.df <- unique(seq.df)
        pydeg_sub <- merge(pydeg_sub,
                           dplyr::select(seq.df,all_of(c("exon_rank",
                                                         "TxP_start",
                                                         "TxP_end",
                                                         "TxP_end",
                                                         "Pseq",
                                                         "indx"))),
                           by = "indx",all.x = TRUE)[,-1]
        
                                        # Write Fasta
        attach(pydeg_sub)
        gr.tmp <- GRanges(seqnames = tx_name,
                          ranges = IRanges(start=TxP_start,
                                           end = TxP_end))
        names(gr.tmp)<-paste(tx_name,
                             comparison,
                             paste0("Peak_coors_G:",peak_start,"-",peak_stop),
                             paste0("Pregion_coors_Tx:",TxP_start,"-",TxP_end))
        mcols(gr.tmp)$sequence <- Pseq

        write2FASTA(gr.tmp,
                    file = file.path(peakSeq_dir,fasta.f),
                    width = 50)
        detach(pydeg_sub)
    }
    
}
