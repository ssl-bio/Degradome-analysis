## 03-Drawing_Dplots.R

## Description: Draws decay plots for all indentified peaks/genes. Two types of plots are drawn, one at the transcript level and another focusing on the peak area (40nt width, user-defined). Checks if the plots exist to avoid overwritting.


##Libraries
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(mgsub))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(RVenn))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PostPyDeg))
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

                                        #Degradome BigWig
dg_bigwig_all <- list()
for (i.sample in sample_list){
                                        #Load bigWig
    bigwig_f <- import(file.path(bigwig_dir,
                                 paste0(i.sample,".mapped_genome_f_CPM.bw")))
    strand(bigwig_f) <- "+"

    bigwig_r <- import(file.path(bigwig_dir,
                                 paste0(i.sample,".mapped_genome_r_CPM.bw")))
    strand(bigwig_r) <- "-"

    dg_bigwig_all[[i.sample]] <- c(bigwig_f,bigwig_r)
}
names(dg_bigwig_all) <- names(sample_list)

                                        #Define mart
mart <- useMart(dataset = "athaliana_eg_gene", biomart = "plants_mart", 
                host = "plants.ensembl.org")

At_genome_seq <- readDNAStringSet(env$At_genome)
##--------------------------------------------------

##Transcript
txdb <- loadDb(file.path(env["supp_data_dir"],"R/txdb_object"))

                                        # Transcript range
tr_range <- GenomicFeatures::transcripts(txdb)
##--------------------------------------------------
##Allow for arbitrary chromosome identifiers
## options(ucscChromosomeNames=FALSE)

for(i.MF in rev(MF_list)) {
    for(i.conf in rev(conf_list)) {
        i.conf_f <- gsub("\\.","_",i.conf)
        ##Only select certain combinations of i.MF and i.conf
        if (i.MF == 4 && i.conf == 0.95 ||
            i.MF == 3 && i.conf == 0.99) {
                                        #Check if file exists
            pydeg_input_f <- file.path(pydeg_pooled_dir,paste("Pooled",
                                                              i.conf_f, "4", i.MF, sep="_"))
            if (file.exists(pydeg_input_f)) {
                
                                        #Read pydeg input file
                pydeg_all <- fread(pydeg_input_f)
                pydeg_all$comparison <- as.factor(pydeg_all$comparison)
                ##Loop over all comparisons
                pydeg_all_list <- list()
                for(i.comparison in levels(pydeg_all$comparison)) {
                                        #Subset comparison
                    pydeg_sub <- pydeg_all[comparison==i.comparison,]

                    if(!is.null(pydeg_sub) && nrow(pydeg_sub)>0) {
                                        #Get comparison between samples
                        i.samples <- unlist(strsplit(gsub("_and_","-",i.comparison),
                                                     split="-"))
                        i.samples2 <- NULL
                        for (i in i.samples) {
                            i.sel <-  sample_list %in% i
                            names(i) <- names(sample_list)[i.sel]
                            i.samples2 <- c(i.samples2,i)
                        }
                        
                        ##control samples
                        i.pairs.ctrl <- unique(names(i.samples2[c(2,4)]))
                        
                        ##test samples
                        i.pairs.test <- unique(names(i.samples2[c(1,3)]))
                        ##-------------------------
                        plot_subdir_gene <- file.path(pydeg_dplot_dir,
                                                      paste0("Gene_", i.comparison))
                        dir.create(plot_subdir_gene)
                        plot_subdir_peak <- file.path(pydeg_dplot_dir,
                                                      paste0("Peak_", i.comparison))
                        dir.create(plot_subdir_peak)
                        ##-----------------------------------
                        ##Subset bam files
                        tmp <- sapply(i.samples,
                                      function(x) {
                                          as.numeric(grepl(x, names(dg_bam_list)))
                                      })
                        
                        sel_bam <- as.logical(apply(tmp,1,sum))
                        
                        bigwigs <- dg_bigwig_all[sel_bam]
                        ##-----------------------------------
                        ##names for labeling
                        sample_list_degradome <-  names(bigwigs)
                        
                                        #match samples to colors
                        names(sample_color) <- sample_list_degradome
                        ##-----------------------------------
                        
                        ##reorder columns
                        ## pydeg_sub <- dplyr::select(pydeg_sub,all_of(i.cols.sort))

                        ## Loop over all combinations of categories
                        i.cat1 <- factor(pydeg_sub$category_1,levels=c(1:4,0))
                        i.cat2 <- as.factor(pydeg_sub$category_2)
                        cat_df <- expand.grid(levels(i.cat1),
                                              levels(i.cat2))
                        
                        for (j.row in seq_len(nrow(cat_df))) {
                            pydeg_plot <- pydeg_sub[pydeg_sub$category_1 == cat_df[j.row,"Var1"],]
                            pydeg_plot <- pydeg_plot[pydeg_plot$category_2 == cat_df[j.row,"Var2"],]

                            if (!is.null(pydeg_plot) && nrow(pydeg_plot)>0) {
                                cat_label <- paste(cat_df[j.row,"Var1"],
                                                   cat_df[j.row,"Var2"],sep = "-")

                                ## Focus on the combination Cat1=1 and Cat2=A
                                if(unique(pydeg_plot$category_1)=="1" &&
                                   unique(pydeg_plot$category_2)=="A") {
                                    Nplots <- length(pydeg_plot$category_1)
                                } else {
                                    Nplots <- 5
                                }
                                
                                sel.na <- is.na(pydeg_plot$gene_region_start)
                                if(sum(sel.na, na.rm=TRUE) > 0) {
                                    pydeg_plot.na <- pydeg_plot[sel.na,]
                                    for (i in seq_len(nrow(pydeg_plot.na))) {
                                        p_start <- pydeg_plot.na[i,"peak_start"]
                                        p_stop <- pydeg_plot.na[i,"peak_stop"]
                                        
                                        max_peak_distace <- max(p_stop) - min(p_start)
                                        plot_offset <- max(40, round((min_plot_width - 
                                                                      max_peak_distace)/2))
                                        pydeg_plot.na$gene_region_start <- max(0, min(p_start) - plot_offset)
                                        pydeg_plot.na$gene_region_end <- max(p_stop) + plot_offset
                                    }

                                    pydeg_plot.p1 <- pydeg_plot[!sel.na,]
                                    pydeg_plot <- rbind(pydeg_plot.p1,pydeg_plot.na)
                                }
                                setDF(pydeg_plot)

                                ##Sort dataframe by txRatio
                                pydeg_plot <- pydeg_plot[with(pydeg_plot,
                                                              order(category_1,
                                                                    category_2,
                                                                    -ratioPTx)),]
                                
                                if (nrow(pydeg_plot) < Nplots) {
                                    np <- nrow(pydeg_plot)
                                } else {
                                    np <- Nplots
                                }

                                ##subset data frame for plotting (should reduce RAM usage)
                                df_plot <- data.frame(pydeg_plot[1:np,])
                                export_pkgs <- c("PostPyDeg", "ggplot2",
                                                 "GenomicRanges", "Cairo",
                                                 "biomaRt", "Gviz")

                                ##Ploting transcript region
                                cl <- makeCluster(as.numeric(env$core))
                                registerDoParallel(cl)
                                foreach(#
                                    i.row = 1:np,
                                    .packages = export_pkgs,
                                    .errorhandling = "remove") %dopar% {
                                        tx.id <- gsub("\\.","_",df_plot[i.row, "tx_name"])
                                        tryCatch({
                                            my.plot <- file.path(#
                                                plot_subdir_gene,
                                                paste0(paste(sprintf("%02d", i.row),"Gene",cat_label,tx.id,
                                                             i.conf_f, "4", i.MF, sep = "_"),
                                                       ".pdf"))
                                            if (!file.exists(my.plot)) {
                                                CairoPDF(my.plot,
                                                         width = 12,
                                                         height = 8)
                                                drawDplot(#
                                                    p_chr = df_plot[i.row, "chr"],
                                                    tx_start = df_plot[i.row, "gene_region_start"],
                                                    tx_end = df_plot[i.row, "gene_region_end"],
                                                    peak_start = df_plot[i.row, "peak_start"], 
                                                    peak_end = df_plot[i.row, "peak_stop"],
                                                    tx_ID = df_plot[i.row, "tx_name"],
                                                    sample_list_plot = sample_list_degradome,
                                                    p_test = i.pairs.test,
                                                    p_ctrl = i.pairs.ctrl,
                                                    p_strand = df_plot[i.row, "strand"],
                                                    ylim = "strand", 
                                                    p_bigwigs = bigwigs,
                                                    p_mart = mart,
                                                    p_colors = sample_color, 
                                                    p_col_hl = col_hl,
                                                    ticksn = 5,
                                                    i_factor=0.15,
                                                    p_width=NULL,
                                                    p_size = 1.8,
                                                    ucscnames=FALSE,
                                                    plot_peak=FALSE,
                                                    p_font=i.font,
                                                    p_base_col=base_col)
                                                dev.off()
                                            }#If there is no plot
                                        }
                                      , error = function(e) {
                                          message(#
                                              paste0(#
                                                  "plotBigWig failed to complete plotting for: chr=", 
                                                  df_plot[i.row, "chr"],
                                                  " gene_region_start=", gene_region_start, 
                                                  " gene_region_end=", gene_region_end,
                                                  " peak_start=", peak_start, 
                                                  " peak_end=", peak_end))
                                          message(e)
                                      })#Try catch
                                    }#[Plot] For each
                                stopCluster(cl)
                                gc()

                                ##Ploting peak region
                                cl <- makeCluster(as.numeric(env$core))
                                registerDoParallel(cl)
                                foreach(#
                                    i.row = 1:np,
                                    .packages = export_pkgs,
                                    .errorhandling = "remove") %dopar% {
                                        tx.id <- gsub("\\.","_",df_plot[i.row, "tx_name"]) 
                                            tryCatch({
                                                my.plot <- file.path(#
                                                    plot_subdir_peak,
                                                    paste0(paste(sprintf("%02d", i.row),"Peak",cat_label,tx.id,
                                                                 i.conf_f, "4", i.MF, sep = "_"),
                                                           ".pdf"))
                                                if (!file.exists(my.plot)) {
                                                    CairoPDF(my.plot,
                                                             width = 12,
                                                             height = 8)
                                                    drawDplot(#
                                                        p_chr = df_plot[i.row, "chr"],
                                                        tx_start = df_plot[i.row, "peak_start"],
                                                        tx_end = df_plot[i.row, "peak_stop"],
                                                        peak_start = df_plot[i.row, "peak_start"], 
                                                        peak_end = df_plot[i.row, "peak_stop"],
                                                        tx_ID = df_plot[i.row, "tx_name"],
                                                        sample_list_plot = sample_list_degradome,
                                                        p_test = i.pairs.test,
                                                        p_ctrl = i.pairs.ctrl,
                                                        p_strand = df_plot[i.row, "strand"],
                                                        ylim = "strand", 
                                                        p_bigwigs = bigwigs,
                                                        p_mart = mart,
                                                        p_colors = sample_color, 
                                                        p_col_hl = col_hl,
                                                        ticksn = 5,
                                                        i_factor=NULL,
                                                        p_width=40,
                                                        p_size = 1.8,
                                                        ucscnames=FALSE,
                                                        plot_peak=TRUE,
                                                        p_font=i.font,
                                                        p_base_col=base_col,
                                                        ref_genome_seq=At_genome_seq,
                                                        ref_genome="TAIR10")
                                                    dev.off()
                                                }#If there is no plot
                                            }
                                          , error = function(e) {
                                              message(#
                                                  paste0(#
                                                      "plotBigWig failed to complete plotting for: chr=", 
                                                      df_plot[i.row, "chr"],
                                                      " gene_region_start=", gene_region_start, 
                                                      " gene_region_end=", gene_region_end,
                                                      " peak_start=", peak_start, 
                                                      " peak_end=", peak_end))
                                              message(e)
                                          })#Try catch
                                        }#[Plot] For each
                                stopCluster(cl)
                                gc()
                            }
                        }
                    }# if isn't null pydeg_sub
                }# Loop over comparions
            }# Test if file with peak candidates exists
        }# Test only a subset of MF and conf levels
    }# Loop over different levels of significance (conf)
}# Loop over different levels of multiplicative factor (MF)

