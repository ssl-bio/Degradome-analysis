##Description draws decay plots for all indentified peaks/genes. Checks if the plot exists to avoid overwritting.
##Panel labels are modified to allow for italics and greek letters but it requires Cairo library for writting the pdf file
##Pending: The function used for drawing is Dplot which draws (barplot) each sample in a different panel perhaps it can be changed to use DplotOL which overlays the replicates (dot plots) to occupy less vertical space

##Libraries
library(here)

library(dplyr)
library(stringr)
## library(reshape2)
library(mgsub)
library(magrittr)
library(data.table)
library(DT)
library(rtracklayer)
library(purrr)
library(RVenn)
library(ggplot2)
library(biomaRt)
library(Gviz)
library(Cairo)
library(GenomicFeatures)
library(Biostrings)
library(RColorBrewer)
library(doParallel)
library(cid7degradomeR)
                                        #library(BSgenome.Athaliana.TAIR.TAIR9)
                                        #library(EnsDb.Athaliana.v51)

##Load environmetnal settings and set paths
if(!exists("env")) {
    library(here)
    here::i_am("Degradome-Oliver-2022-R/00-Initialization.R")
    source(here("Degradome-Oliver-2022-R/00-Initialization.R"))
}

##--------------------------------------------------
                                        #Degradome BigWig
dg_bigwig_all <- list()
for (i.sample in sample_list){
                                        #Load bigWig
    bigwig_f <- import(file.path(bigwig_dir,paste0(i.sample,".mapped_genome_f_CPM.bw")))
    strand(bigwig_f) <- "+"

    bigwig_r <- import(file.path(bigwig_dir,paste0(i.sample,".mapped_genome_r_CPM.bw")))
    strand(bigwig_r) <- "-"

    dg_bigwig_all[[i.sample]] <- c(bigwig_f,bigwig_r)
}
names(dg_bigwig_all) <- names(sample_list)

                                        #Plot settings
col_cond <- c("#FC4E2A", "#E31A1C", "#2c84a8ff", "#214b7fff")

                                        #Color for highlight in profile plot
col_hl <- "#e5ff00bf" #f2ff00cc#e5ff00d9
                                        #"#00FF0099"#ffef00ff#ffef00bf#ffef0080
base_col <- c( A="#33A02C",C="#1F78B4",U="#E31A1C", G="#FF7F00")
                                        #font
i.font <- ps.options()$family
                                        #minimum width of plot
min_plot_width <- 1000

                                        #Define mart
mart <- useMart(dataset = "athaliana_eg_gene", biomart = "plants_mart", 
                host = "plants.ensembl.org")

At_genome <- readDNAStringSet(env$At_genome)
##--------------------------------------------------

##Transcript
txdb <- loadDb(file.path(env["supp_data_dir"],"txdb_object"))

                                        # Transcript range
tr_range <- GenomicFeatures::transcripts(txdb)
##--------------------------------------------------
i.tx <- c()

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

                sel.tx <- pydeg_all$tx_name %in% i.tx
                pydeg_all <- pydeg_all[sel.tx,]
                
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
                        plot_subdir_gene <- file.path(pydeg_dplot_dir,paste0("Gene_",i.comparison))
                        dir.create(plot_subdir_gene)
                        plot_subdir_peak <- file.path(pydeg_dplot_dir,paste0("Peak_",i.comparison))
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
                        names(col_cond) <- sample_list_degradome
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
                                export_pkgs <- c("cid7degradomeR", "ggplot2",
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
                                                    p_colors = col_cond, 
                                                    p_col_hl = col_hl,
                                                    ticksn = 5,
                                                    i_factor=0.15,
                                                    p_width=NULL,
                                                    p_size = 1.8,
                                                    ucscnames=FALSE,
                                                    plot_peak=FALSE,
                                                    p_font=i.font,
                                                    p_base_col=base_col,
                                                    ref_genome=At_genome)
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
                                                        p_colors = col_cond, 
                                                        p_col_hl = col_hl,
                                                        ticksn = 5,
                                                        i_factor=NULL,
                                                        p_width=40,
                                                        p_size = 1.8,
                                                        ucscnames=FALSE,
                                                        plot_peak=TRUE,
                                                        p_font=i.font,
                                                        p_base_col=base_col,
                                                        ref_genome=At_genome)
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

