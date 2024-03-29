## 03-Drawing_Dplots.R

## Description: Draws decay plots for all indentified peaks/genes. Two types of plots are drawn, one at the transcript level and another focusing on the peak area (40nt width, user-defined). Checks if the plots exist to avoid overwritting.


## Libraries
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mgsub))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(RVenn))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(PostPyDeg))
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

## Minimal set variables
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

## Degradome BigWig
dg_bigwig_all <- list()
for (i.sample in sample_list){
    ## Load bigWig
    bigwig_f <- import(file.path(bigwig_dir,
                                 paste0(i.sample, "_G_f_DESeq.bw")))
    strand(bigwig_f) <- "+"

    bigwig_r <- import(file.path(bigwig_dir,
                                 paste0(i.sample, "_G_r_DESeq.bw")))
    strand(bigwig_r) <- "-"

    dg_bigwig_all[[i.sample]] <- c(bigwig_f, bigwig_r)
}
names(dg_bigwig_all) <- names(sample_list)

## Define mart
mart <- useMart(dataset = "athaliana_eg_gene", biomart = "plants_mart",
                host = "https://plants.ensembl.org")

At_genome_seq <- readDNAStringSet(env$At_genome)
##--------------------------------------------------

## Transcript
txdb <- loadDb(file.path(env["supp_data_dir"], "R/txdb_object"))

## Transcript range
tr_range <- GenomicFeatures::transcripts(txdb)

## Default number of plots per combination of categories
Nplots_default <- 50
Nplots_max <- 100
##--------------------------------------------------

for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    i.conf_f <- gsub("\\.", "_", i.conf)
    input_file <- file.path(pydeg_pooled_dir,
                            paste("Pooled", i.conf_f, "4", i.MF, sep = "_"))
    if (!file.exists(input_file)) {
        next
    }

    ## Read pydeg input file
    pydeg_all <- fread(input_file)
    pydeg_all$comparison <- as.factor(pydeg_all$comparison)

    ##Loop over all comparisons
    pydeg_all_list <- list()
    for (i.comparison in levels(pydeg_all$comparison)) {
        ##Subset comparison
        pydeg_sub <- pydeg_all[comparison == i.comparison,]

        if (is.null(pydeg_sub) || nrow(pydeg_sub)==0) {
            next
        }
        i.samples <- unlist(strsplit(gsub("_and_", "-", i.comparison),
                                     split = "-"))
        i.samples2 <- NULL
        for (i in i.samples) {
            i.sel <-  sample_list %in% i
            names(i) <- names(sample_list)[i.sel]
            i.samples2 <- c(i.samples2, i)
        }

        ##control samples
        i.pairs.ctrl <- unique(names(i.samples2[c(2, 4)]))

        ##test samples
        i.pairs.test <- unique(names(i.samples2[c(1, 3)]))
        ##-------------------------
        plot_subdir_gene <- file.path(pydeg_dplot_dir,
                                      paste0("Gene_", i.comparison))
        dir.create(plot_subdir_gene, recursive=TRUE)
        plot_subdir_peak <- file.path(pydeg_dplot_dir,
                                      paste0("Peak_", i.comparison))
        dir.create(plot_subdir_peak, recursive=TRUE)
        ##-----------------------------------
        ##Subset bam files
        tmp <- sapply(i.samples,
                      function(x) {
                          as.numeric(grepl(x, names(dg_bam_list)))
                      })

        sel_bam <- as.logical(apply(tmp, 1,sum))

        bigwigs <- dg_bigwig_all[sel_bam]
        ##-----------------------------------
        ##names for labeling
        sample_list_degradome <-  names(bigwigs)

        ## match samples to colors
        names(sample_color) <- sample_list_degradome
        ##-----------------------------------
        ## Loop over all combinations of categories
        i.cat1 <- factor(pydeg_sub$category_1, levels = c(1:4, 0))
        i.cat2 <- as.factor(pydeg_sub$category_2)
        cat_df <- expand.grid(levels(i.cat1),
                              levels(i.cat2))

        for (j.row in seq_len(nrow(cat_df))) {
            pydeg_plot <- pydeg_sub[pydeg_sub$category_1 == cat_df[j.row, "Var1"], ]
            pydeg_plot <- pydeg_plot[pydeg_plot$category_2 == cat_df[j.row, "Var2"], ]

            if (is.null(pydeg_plot) || nrow(pydeg_plot) == 0) {
                next
            }
            cat_label <- paste(cat_df[j.row, "Var1"],
                               cat_df[j.row, "Var2"], sep = "-")

            ## Focus on the combination Cat1 = 1 and Cat2 = A
            if (unique(pydeg_plot$category_1) == "1" &&
                unique(pydeg_plot$category_2) == "A") {
                if (length(pydeg_plot$category_1) > Nplots_max) {
                    Nplots <- Nplots_max
                } else {
                    Nplots <- length(pydeg_plot$category_1)
                }
            } else {
                if (length(pydeg_plot$category_1) < Nplots_default) {
                    Nplots <- length(pydeg_plot$category_1)
                } else {
                    Nplots <- Nplots_default
                }
            }

            sel.na <- is.na(pydeg_plot$gene_region_start)
            if (sum(sel.na, na.rm = TRUE) > 0) {
                pydeg_plot.na <- pydeg_plot[sel.na, ]
                for (i in seq_len(nrow(pydeg_plot.na))) {
                    p_start <- pydeg_plot.na[i, "peak_start"]
                    p_stop <- pydeg_plot.na[i, "peak_stop"]

                    max_peak_distace <- max(p_stop) - min(p_start)
                    plot_offset <- max(40, round((min_plot_width -
                                                  max_peak_distace) / 2))
                    ## pydeg_plot.na$gene_region_start <- max(0, min(p_start) - plot_offset)
                    ## pydeg_plot.na$gene_region_end <- max(p_stop) + plot_offset
                    pydeg_plot.na <- pydeg_plot.na |>
                        mutate(gene_region_start = max(0, min(p_start) - plot_offset)) |>
                        mutate(gene_region_end <- max(p_stop) + plot_offset)
                }

                pydeg_plot.p1 <- pydeg_plot[!sel.na, ]
                pydeg_plot <- rbind(pydeg_plot.p1, pydeg_plot.na)
            }
            setDF(pydeg_plot)

            ##Sort dataframe by txRatio
            pydeg_plot <- pydeg_plot |>
                arrange(category_1, category_2, desc(ratioPTx))

            if (nrow(pydeg_plot) < Nplots) {
                np <- nrow(pydeg_plot)
            } else {
                np <- Nplots
            }

            ##subset data frame for plotting (should reduce RAM usage)
            df_plot <- data.frame(pydeg_plot[1:np, ])
            export_pkgs <- c("PostPyDeg", "ggplot2",
                             "GenomicRanges", "Cairo",
                             "biomaRt", "Gviz")

            ## Ploting transcript region
            my.gene.plots <- file.path(
                plot_subdir_gene,
                paste0(
                    paste(
                        sprintf("%02d",
                                seq(1,
                                    nrow(df_plot))),
                        "Gene",
                        cat_label,
                        gsub("\\.",
                             "_", df_plot$tx_name),
                        i.conf_f, "4",
                        i.MF, sep = "_"),
                    ".pdf")
            )

            gene.sel <- !file.exists(my.gene.plots)

            if (sum(gene.sel) == 0) {
                next
            }
            df_plot$file <- my.gene.plots
            df_plot_gene <- df_plot[gene.sel, ]
            cl <- makeCluster(as.numeric(env$core))
            registerDoParallel(cl)
            foreach(#
                i.row = seq_len(nrow(df_plot_gene)),
                .packages = export_pkgs,
                .errorhandling = "remove") %dopar% {
                    tx.id <- gsub("\\.", "_", df_plot_gene[i.row,  "tx_name"])
                    ## Ploting transcript region
                    tryCatch({
                        cairo_pdf(file = df_plot_gene$file[i.row],
                                  width = 9,
                                  height = 6)
                        drawDplot(#
                            p_chr = df_plot_gene[i.row, "chr"],
                            tx_start = df_plot_gene[i.row, "gene_region_start"],
                            tx_end = df_plot_gene[i.row, "gene_region_end"],
                            peak_start = df_plot_gene[i.row, "peak_start"],
                            peak_end = df_plot_gene[i.row, "peak_stop"],
                            tx_ID = df_plot_gene[i.row, "tx_name"],
                            sample_list_plot = sample_list_degradome,
                            p_test = i.pairs.test,
                            p_ctrl = i.pairs.ctrl,
                            p_strand = df_plot_gene[i.row, "strand"],
                            ylim = "strand",
                            p_bigwigs = bigwigs,
                            p_mart = mart,
                            p_colors = sample_color,
                            p_col_hl = col_hl,
                            ticksn = 5,
                            i_factor = 0.15,
                            p_width = NULL,
                            p_size = 1.8,
                            ucscnames = FALSE,
                            plot_peak = FALSE,
                            p_font = i.font,
                            p_base_col = base_col)
                        dev.off()
                    }
                  , error = function(e) {
                      message(#
                          paste0(#
                              "plotBigWig failed to complete plotting for: chr = ",
                              df_plot_gene[i.row, "chr"],
                              " gene_region_start = ", gene_region_start,
                              " gene_region_end = ", gene_region_end,
                              " peak_start = ", peak_start,
                              " peak_end = ", peak_end))
                      message(e)
                  })#Try catch
                }#[Plot] For each
            stopCluster(cl)
            gc()

            ## Ploting peak region
            my.peak.plots <- file.path(plot_subdir_peak,
                                       paste0(
                                           paste(
                                               sprintf("%02d",
                                                       seq(1,
                                                           nrow(df_plot))),
                                               "Peak",
                                               cat_label,
                                               gsub("\\.", "_",
                                                    df_plot$tx_name),
                                               i.conf_f, "4", i.MF, sep = "_"),
                                           ".pdf"))

            peak.sel <- !file.exists(my.peak.plots)

            if (sum(peak.sel) == 0) {
                next
            }
            df_plot$file <- my.peak.plots
            df_plot_peak <- df_plot[peak.sel, ]
            cl <- makeCluster(as.numeric(env$core))
            registerDoParallel(cl)
            foreach(#
                i.row = seq_len(nrow(df_plot_peak)),
                .packages = export_pkgs,
                .errorhandling = "remove") %dopar% {
                    tx.id <- gsub("\\.", "_",
                                  df_plot_peak[i.row,
                                               "tx_name"])
                    ## Ploting peak region
                    tryCatch({
                        cairo_pdf(file = df_plot_peak$file[i.row],
                                  width = 9,
                                  height = 6)
                        drawDplot(#
                            p_chr = df_plot_peak[i.row, "chr"],
                            tx_start = df_plot_peak[i.row, "peak_start"],
                            tx_end = df_plot_peak[i.row, "peak_stop"],
                            peak_start = df_plot_peak[i.row, "peak_start"],
                            peak_end = df_plot_peak[i.row, "peak_stop"],
                            tx_ID = df_plot_peak[i.row, "tx_name"],
                            sample_list_plot = sample_list_degradome,
                            p_test = i.pairs.test,
                            p_ctrl = i.pairs.ctrl,
                            p_strand = df_plot_peak[i.row, "strand"],
                            ylim = "strand",
                            p_bigwigs = bigwigs,
                            p_mart = mart,
                            p_colors = sample_color,
                            p_col_hl = col_hl,
                            ticksn = 5,
                            i_factor = NULL,
                            p_width = 40,
                            p_size = 1.8,
                            ucscnames = FALSE,
                            plot_peak = TRUE,
                            p_font = i.font,
                            p_base_col = base_col,
                            ref_genome_seq = At_genome_seq,
                            ref_genome = "TAIR10")
                        dev.off()
                    }
                  , error = function(e) {
                      message(#
                          paste0(#
                              "plotBigWig failed to complete plotting for: chr = ",
                              df_plot_peak[i.row, "chr"],
                              " gene_region_start = ",
                              gene_region_start,
                              " gene_region_end = ",
                              gene_region_end,
                              " peak_start = ", peak_start,
                              " peak_end = ", peak_end))
                      message(e)
                  })#Try catch
                }#[Plot] For each
            stopCluster(cl)
            gc()
        }
        ## }# if isn't null pydeg_sub
    }# Loop over comparions
    ## }# Test if file with peak candidates exists
}
