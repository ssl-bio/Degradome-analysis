## 04-getPeakSeq.R

## Description: From candidate peaks classified into groups 1-A, it produces a fasta file with sequences around the peak (44nt window). Transcript name, comparison where peak was found and genomic and transcrip coordinates are included in the header of each fasta entry

## Output dir: Supporting_data/Peak_sequences
## Output files:
## - PeakRegioncDNA_<conf>_<win>_<mf>.fa

## Example:
## - PeakRegioncDNA_0_95_4_2.fa

## Libraries
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

## Parse arguments
option_list <- list(
    make_option(c("-d", "--wd"), type = "character", default = NULL,
                help = "Working directory", metavar = "character"),
    make_option(c("-b", "--base"), type = "character", default = NULL,
                help = "Project basename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

setwd(opt$wd)

## Import system variables
ivarF <- "Initialization_variables.RData"
ivarM <- paste0("minimal_variables_", opt$base, ".RData")
error_msg <- " was not found.\nPlease run 'Scripts/R/00-Initialization.R' providing the root directory and the project base name.\nAlternatively run '02-Degradome.sh' providing the project base name.\n"

## Minimal set of variables
min_variables <- file.path("Env_variables", ivarM)
if (file.exists(min_variables)) {
    load(min_variables)
    ivars <- file.path(supp_data_dir, "R/Initialization_variables.RData")
    if (file.exists(ivars)) {
        load(ivars)
    } else {
        cat(ivarF, error_msg)
        break
    }
} else {
    cat(ivarM, error_msg)
    break
}
##-------------------------------------------------------------------------

EDB <- EnsDb(DB)

## Create a dataframe with Tx sequences
i.fasta <- readDNAStringSet(env[["At_transcript"]])
names(i.fasta) <-  vapply(strsplit(names(i.fasta), " "),
                          `[`, 1, FUN.VALUE = character(1))
tx_sequences <- data.frame(tx_id=names(i.fasta),
                           tx_width=width(i.fasta),
                           seq=as.character(i.fasta))

## fasta_names <- names(i.fasta)
peak_region <- 44

## Create a dataframe with unique peaks from Cat 1-A 1-B 2-A 2-B
pooled_list <- list()
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    i.conf_f <- gsub("\\.", "_", i.conf)

    ## Define input and output files
    input_file <- file.path(pydeg_pooled_dir, paste("Pooled",
                                                    i.conf_f, "4",
                                                    i.MF, sep = "_"))
    pydeg_all <- fread(input_file)
    sel <- (pydeg_all$category_1 == "1" |  pydeg_all$category_1 == "2") &
        (pydeg_all$category_2 == "A" | pydeg_all$category_2 == "B")
    
    if (sum(sel) == 0) {
        ## sel <- (pydeg_all$category_1 == "3" |  pydeg_all$category_1 == "4")
        next
    }
    pydeg_sub <- pydeg_all[sel,] |>
        mutate(indx = paste0(tx_name, peak_start,
                             peak_stop)) |>
        dplyr::select(all_of(c("tx_name", "chr",
                               "strand", "peak_start",
                               "peak_stop","indx")))
    
    pooled_list[[i]] <- pydeg_sub
    
}
pooled_dt <- do.call(rbind, pooled_list)
if (!is.null(pooled_dt) && nrow(pooled_dt) > 0) {
    isel <- !duplicated(pooled_dt)
    pooled_dt <- pooled_dt[isel, ]

    ## merge dataframes
    gnm <- GRanges(pooled_dt$chr,
                   IRanges(start = pooled_dt$peak_start,
                           end = pooled_dt$peak_stop),
                   strand = pooled_dt$strand)
    gnm_tx <- genomeToTranscript(gnm, EDB)
    names(gnm_tx) <- pooled_dt$indx
    gnm_tx_df <- as.data.frame(gnm_tx)
    gnm_tx_df <- merge(gnm_tx_df, tx_sequences, by = "tx_id", all.x = TRUE)
    setnames(gnm_tx_df,"group_name","indx")

    ## Loop over pydeg settings
    for (i in seq_along(MF_list)) {
        i.MF <- MF_list[i]
        i.conf <- conf_list[i]
        i.conf_f <- gsub("\\.", "_", i.conf)

        ## Define input and output files
        input_file <- file.path(pydeg_pooled_dir, paste("Pooled",
                                                        i.conf_f, "4",
                                                        i.MF, sep = "_"))
        out_file <- file.path(peakSeq_dir,
                              paste0(paste("PeakRegioncDNA",
                                           i.conf_f, "4", i.MF,
                                           sep = "_"), ".fa"))

        if (file.exists(out_file) || !file.exists(input_file)) {
            next
        }

        pydeg_all <- fread(input_file)
        sel <- (pydeg_all$category_1 == "1" |  pydeg_all$category_1 == "2") &
            (pydeg_all$category_2 == "A" | pydeg_all$category_2 == "B")

        if (sum(sel) == 0) {
            ## sel <- (pydeg_all$category_1 == "3" |  pydeg_all$category_1 == "4")
            next
        }

        pydeg_sub <- pydeg_all[sel, ] |>
            mutate(indx = paste0(tx_name, peak_start, peak_stop))

        ## Merge to get tx coordinates and sequences
        pydeg_sub_tx  <- merge(pydeg_sub,
                               dplyr::select(gnm_tx_df,
                                             all_of(c("tx_id",
                                                      "tx_width",
                                                      "start",
                                                      "end",
                                                      "exon_rank",
                                                      "seq",
                                                      "indx"))),
                               by = "indx", all.x = TRUE)[, -1]
        isel <- pydeg_sub_tx$tx_name == pydeg_sub_tx$tx_id
        pydeg_sub_tx <- pydeg_sub_tx[isel, ]

        ## Get transcript coordinates to extract sequences around peaks
        width_modulo <- (pydeg_sub_tx$peak_width+1) %% 2
        seq_offset <- (peak_region - (pydeg_sub_tx$peak_width + 1) - width_modulo) / 2
        
        pydeg_sub_tx <- pydeg_sub_tx |>
            mutate(i.start = start - seq_offset) |>
            mutate(i.end = end + seq_offset + width_modulo) |>
            mutate(i.start = if_else(i.start < 1,
                                     1, i.start)) |>
            mutate(i.end = if_else(i.end > tx_width.y,
                                   tx_width.y, i.end)) |>
            mutate(peak_seq = subseq(seq,
                                     start = i.start,
                                     end = i.end))
        
        ## Write Fasta file
        ## attach(pydeg_sub_tx)
        gr.tmp <- GRanges(seqnames = pydeg_sub_tx$tx_name,
                          ranges = IRanges(start = pydeg_sub_tx$i.start,
                                           end = pydeg_sub_tx$i.end))
        names(gr.tmp) <- paste(pydeg_sub_tx$tx_name,
                               pydeg_sub_tx$comparison,
                               paste0("Peak_coors_G:",
                                      pydeg_sub_tx$peak_start, "-",
                                      pydeg_sub_tx$peak_stop),
                               paste0("Pregion_coors_Tx:",
                                      pydeg_sub_tx$i.start, "-",
                                      pydeg_sub_tx$i.end))
        mcols(gr.tmp)$sequence <- pydeg_sub_tx$peak_seq

        write2FASTA(gr.tmp,
                    file = out_file,
                    width = 50)
    }
}
