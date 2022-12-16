## 00-Initialization.R

## Description: Script to set up most of the variables to be used in the rest of the analysis,  creates the directories to store the output files and also produces auxiliary files to speed up the annotation of peaks.

## Output dirs created:
## [PyDegradome output with annotations]
## <output_root_dir>/output_02/01-PyDegradome_processed

## [PyDegradome_processed output classified and combined]
## <output_root_dir>/output_02/02-PyDegradome_pooled

## [Summary counts of all classification steps]
## <output_root_dir>/output_02/03-Report/Summary

## [Tables with mapping reads to speed up process]
## <output_root_dir>/Supporting_data/Degradome_reads

## [Tables with coordinates of non-peak regions]
## <output_root_dir>/Supporting_data/PyDegradome_NonPeaks_coordinates

## [Tables with scaled mapping reads (test peaks)]
## <output_root_dir>/Supporting_data/PyDegradome_maxPeak_Scaled

## [Tables with scaled mapping reads (ctrl transcript)]
## <output_root_dir>/Supporting_data/PyDegradome_maxTx_Scaled

## [Fasta files with genomic sequences around peaks]
## <output_root_dir>/Supporting_data/Peak_sequences

## [Auxilary files to speed up running of scripts]   
## <output_root_dir>/Supporting_data/R

## Auxiliary files created:
## [Variables defined in this script including those in the Variable definition file]
## - Initialization_variables.RData

## [A Sqlite database used to transform genomic coordinates to transcript coordinates]
## - Arabidopsis_thaliana.TAIR10.55.sqlite

## [Table from plant ensembl to obtain the longest isoform (representative gene)]
## - Transcript_information.txt

## [A transcript database object build from the annotation file]
## - txdb_object

## Libraries
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(filesstrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ensembldb))
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
## Path settings
env.raw <- read.table(paste0("Env_variables/Degradome_",
                             opt$base, ".txt"),
                      sep = "=", quote = "")
env.raw[,2] <- gsub('" "', '"\\, "', env.raw[, 2])
isel <- grepl("\\(", env.raw[, 2])
env.raw1 <- env.raw[isel, ]
env.raw2 <- env.raw[!isel, ]

for (i in seq_len(nrow(env.raw1))) {
    istr <- env.raw1[i, 2]
    assign(env.raw1[i, 1],
           eval(parse(text = paste0("c",
                                  gsub('" "', '"\\, "', istr)))))
}

env <- str_trim(env.raw2$V2)
names(env) <- str_trim(env.raw2$V1)
env <- as.list(env)



## If `core` is not defined in env, use a single core.
if (!is.element("cores", names(env))) {
    env$cores <- 1
} else {
    env$cores <- as.numeric(env$cores)
}

## Input path
                                        # output_1 list directories
idirs <- list.dirs(env$output_dirB,
                   recursive  = FALSE,
                   full.names = FALSE)

                                        # assign selected directories
input_df <- data.frame(ivars = c("bam_dir", "bam_tr_dir",
                               "bigwig_dir", "bigwig_dir", "pydeg_dir",
                               "htseq_dir"),
                       ikeys = c("bam_genomic", "bam_transcript",
                               "bigwig_genomic", "bigwig_chromosome",
                               "PyDegradome",  "htseq_genomic"))

for (i in seq_len(nrow(input_df))) {
    ivar <- input_df[i, 1]
    ikey <- input_df[i, 2]
    isel <- grepl(ikey, idirs)
    if(sum(isel) == 1) {
        assign(ivar,
               eval(file.path(env["output_dirB"], idirs[isel])))
    }
}

## Output path
pydeg_processed_dir <- file.path(env["output_dirR"],
                                 "01-PyDegradome_processed")
pydeg_pooled_dir <- file.path(env["output_dirR"],
                              "02-PyDegradome_pooled")
report_dir <- file.path(env["output_dirR"],
                        "03-Report")
summary_dir <- file.path(report_dir,
                         "Summary")

## Peak reads
pydeg_reads_dir <- file.path(env["supp_data_dir"],
                             "Degradome_reads")

## Coordinates of regions with no peaks
pydeg_np_dir <- file.path(env["supp_data_dir"],
                          "PyDegradome_NonPeaks_coordinates")

## Max peak reads
maxP_dir <- file.path(env["supp_data_dir"],
                      "PyDegradome_maxPeak_Scaled")

## Max transcript reads
maxR_dir <- file.path(env["supp_data_dir"],
                      "PyDegradome_maxTx_Scaled")

## Root dir to store degradome plots
pydeg_dplot_dir <- file.path(report_dir,
                             "Dplots")

## Fasta dir
peakSeq_dir <- file.path(env["supp_data_dir"],
                         "Peak_sequences")

## R various
Rvarious_dir <- file.path(env["supp_data_dir"], "R")

dir.create(env[["output_dirR"]])
dir.create(report_dir)
dir.create(pydeg_pooled_dir)
dir.create(pydeg_processed_dir)
dir.create(env[["supp_data_dir"]])
dir.create(summary_dir)
dir.create(pydeg_reads_dir)
dir.create(pydeg_np_dir)
dir.create(maxR_dir)
dir.create(maxP_dir)
dir.create(pydeg_dplot_dir)
dir.create(peakSeq_dir)
dir.create(Rvarious_dir)
                                        #--------------------------------------------------------------------

## Settings used for PyDegradome analysis
pydeg_settings <- unlist(strsplit(pydeg_script_settings, " "))
MF_list <- as.numeric(pydeg_settings[rep(c(FALSE, TRUE),
                                         length(pydeg_script_settings))])
conf_list <- as.numeric(pydeg_settings[rep(c(TRUE, FALSE),
                                           length(pydeg_script_settings))])
## Sample name list
names(control_samples) <- control_samples_name
names(test_samples) <- test_samples_name

idf <- expand.grid(test=test_samples,nc=control_samples)
if (nrow(idf)==1) {
    idf <- rbind(idf,idf)
}
comp_list <- list()
for (i in seq_along(colnames(idf))) {
    i.not <- seq_along(colnames(idf))[!seq_along(colnames(idf)) %in% i]
    i.comp <- colnames(idf)[i]
    i.comp_other <- colnames(idf)[i.not]
    comp_list_sub <- NULL
    for (j in seq_len(nrow(idf))) {
        x <- paste("t", idf[j, i.comp],
                   "c", idf[j, i.comp_other], sep="_")
        comp_list_sub <- c(comp_list_sub, x)
    }
    comp_list[[i.comp]] <- comp_list_sub
}

tmp <- data.frame(sapply(comp_list, c))


tmp2 <- list(as.data.frame(combn(tmp[["test"]], 2)),
             as.data.frame(combn(tmp[["nc"]], 2)))
comp_pair_list <- list()
for (k in seq_along(tmp2[[1]])) {
    comp_pair_list[[k]] <- list(test=paste(tmp2[[1]][, k],
                                           collapse = "_and_"),
                                nc=paste(tmp2[[2]][,k],
                                         collapse = "_and_"))
}

treatments <- gsub(" \\[[0-9]\\]","",
                   c(names(control_samples),names(test_samples)))

idesign <- data.frame(sample = c(control_samples,
                               test_samples),
                      replicate = c(seq_along(control_samples),
                                  seq_along(test_samples)),
                      treatment = treatments)
rownames(idesign) <- seq_len(nrow(idesign))

sample_list <- c(test_samples, control_samples)

dg_bigwig_list <- structure(file.path(bam_dir, sample_list),
                            names = sample_list)

dg_bam_list <- structure(file.path(bam_dir, paste0(sample_list,
                                                  ".bam")),
                         names = sample_list)

dg_bam_tr_list <- structure(file.path(bam_tr_dir,
                                      paste0(sample_list,
                                             "_uni_sort.bam")),
                            names = sample_list)
##--------------------------------------------------
##Build a transcript info data frame
##Data frame with feature, start, stop, gene name information from biomaRt
tx_info_f <- file.path(Rvarious_dir,"Transcript_information.txt")
if (!file.exists(tx_info_f)) {
    sp_split <- unlist(strsplit(env$sp,"_"))
    isp <- paste0(substr(sp_split[1],1,1),sp_split[2],"_eg_gene")
    mart <- useMart(dataset = isp,
                    biomart = "plants_mart",
                    host = "plants.ensembl.org")

    idescription <- getBM(
        mart = mart,
        attributes = c(
            "ensembl_transcript_id",
            "transcript_start",
            "transcript_end",
            "external_gene_name",
            "entrezgene_description",
            "interpro_short_description")
    )

    icoordinates <- getBM(
        mart = mart,
        attributes = c("ensembl_transcript_id",
                       "cds_length"))

    iannot <- merge(idescription, icoordinates,
                    by = 'ensembl_transcript_id', all.x = TRUE)

    colnames(iannot) <- c("tx_name", "gene_region_start", "gene_region_end",
                          "gene_name", "entrezgene_description",
                          "interpro_short_description", "cds_len")
   
    i.cols <- c("ID", "tx_name", "gene_region_start",
                "gene_region_end", "cds_len", "gene_name", "Description")

    iannot_locus_yes <- iannot[iannot$tx_name!="", ]
    iannot_locus_yes$ID <- gsub("\\.[0-9]", "", iannot_locus_yes$tx_name)

                                        # Create a Description column based on the entrezgene and interpro descriptions
    ilog <- iannot_locus_yes$entrezgene_description == ""
    temp1 <- iannot_locus_yes[ilog, ]
    temp2 <- iannot_locus_yes[!ilog, ]

    temp1$Description <- temp1$interpro_short_description
    temp2$Description <- temp2$entrezgene_description

    temp1 <- dplyr::select(temp1, all_of(i.cols))
    temp2 <- dplyr::select(temp2, all_of(i.cols))
    tx_info <- rbind(temp1, temp2)
   
   
    ## Create an index column to remove duplicated entries
    tx_info$indx <-  paste0(tx_info$tx_name,
                            tx_info$gene_region_start,
                            tx_info$gene_region_end)
   
    sel <- duplicated(tx_info$indx)
    tx_info <- dplyr::select(tx_info[!sel,],all_of(i.cols))


    ## tx_info <- tx_info[with(tx_info,order(tx_name)),]

    tx_info$tx_width <- tx_info$gene_region_end - tx_info$gene_region_start

    ## Identify genes with and without isoforms
    n_occur_id <- data.frame(table(tx_info$ID))
    i.log.dup <- tx_info$ID %in% n_occur_id$Var1[n_occur_id$Freq > 1]
    df_dup <- tx_info[i.log.dup,]
    df_uniq <- tx_info[!i.log.dup,]

                                        #work on genes with isoforms
                                        #Add representative gene model
    dup_ID <- unique(df_dup$ID)
    df_np_list <- list()
    for(k in seq_along(dup_ID)) {
        i.id <- dup_ID[k]
        temp_df <- df_dup[df_dup$ID == i.id, ]
        ## tx_max <- max(temp_df$tx_width)
        cds_max <- max(temp_df$cds_len)
        temp_df$rep_gene <- if_else(temp_df$cds_len == cds_max, 1, 0)

        if (length(which(temp_df$rep_gene == 1)) > 1) {
            ## cds length is the same, select the first isoform
            min_tx <- min(temp_df$tx_name)
            temp_df$rep_gene <- if_else(temp_df$tx_name == min_tx, 1, 0)
        }
        df_np_list[[k]] <- temp_df
    }
    df_np <- do.call(rbind, df_np_list)

                                        #work on genes without isoforms
                                        #Add representative gene model
    df_uniq$rep_gene <- 1

                                        #Merge dataframes and sort
    tx_info <- rbind(df_np, df_uniq)
    tx_info <- tx_info[with(tx_info, order(tx_name)), ]
   
    fwrite(tx_info, tx_info_f, quote = FALSE, sep = "\t")
}
##--------------------------------------------------

## Make a TxDb object from transcript annotations available as a GFF3 or GTF file
txdb_f <- file.path(Rvarious_dir, "txdb_object")
if (!file.exists(txdb_f)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(#
                                 file = env[["ref_gtf"]], format = "gtf")
    saveDb(txdb, txdb_f)
}

## Generate a sqlite database
ensemblDB_file <- file.path(Rvarious_dir,
                            gsub("\\.gtf", ".sqlite",
                                 basename(env[["ref_gtf"]])))
if (!file.exists(ensemblDB_file)) {
    DB <- ensDbFromGtf(gtf = env[["ref_gtf"]], path=Rvarious_dir)
} else {
    DB <- ensemblDB_file
}


## Columns to round up
## cols2round <- c("max_non_peak_ratio", "max_non_peak_ratio_2", "max.peak.scaled", "max.read.tx.scaled", "ratioPTx")

## Link to TAIR
tair.prefix <- "<a  target=_blank href=https://www.arabidopsis.org/servlets/TairObject?type=locus&name="
tair.suffix <- ">TAIR</a>"

## Link to PARE
pare.prefix <- "<a  target=_blank href=https://mpss.danforthcenter.org/web/php/pages/GeneAnalysis.php?SITE=at_pare&featureName="
pare.suffix <- "&model=1>"

## Link to imges
plot.prefix <- "<a target=_blank href="
plot.suffix1 <- ' onclick=\"window.open('
plot.suffix_widthG <- ",'geneplot', 'resizable, width=1060, height=758'); return false;"
plot.suffix_widthP <- ",'peakplot', 'resizable, width=1060, height=758'); return false;"
plot.suffix_widthS <- ",'algwindowA', 'resizable, width=575, height=290'); return false;"
plot.suffix_widths <- ",'algwindowB', 'resizable, width=450, height=220'); return false;"

## List of columns for rounding values
cols2round <- c("Peak:NonPeak (1)", "Peak:NonPeak (2)", "Mean Peak Test (Scaled)", "Mean Max Tx Ctrl (Scaled)", "Peak(T):Max(C)", "DeltaG\nopen")

## Columns for exporting table of candidates
i.cols.sort <- c("cat1_2",
                 "tx_name",
                 "Gene_plot",
                 "Peak_plot",
                 "rep_gene",
                 "feature_type",
                 "gene_name",
                 "ratioPTx",
                 ## "category_1",
                 ## "category_2",
                 "MorePeaks",
                 "Description",
                 "comparison",
                 "max_peak_1",#Metrics
                 "max_np_gene_1",
                 "max_non_peak_ratio_1",
                 "max_peak_2",
                 "max_np_gene_2",
                 "max_non_peak_ratio_2",
                 "shared",
                 "max_peak_scaled",
                 "max_read_tx_scaled",
                 "chr",#Coordinates
                 "strand",
                 "peak_start",
                 "peak_stop",
                 "width",
                 "gene_region_start",
                 "gene_region_end",
                 "feature_start",
                 "feature_end",
                 "cds_len",
                 "tx_width")

i.cols.export <- c("Cat\n1&2",#Main
                   "Transcript",
                   "Gene\nplot",
                   "Peak\nplot",
                   "Is rep.",
                   "Feature",
                   "Gene name",
                   "Peak(T):Max(C)",
                   ## "Cat. 1",
                   ## "Cat. 2",
                   "Peaks 2>",
                   "Description",
                   "Comparison",
                   "Max Peak (1)", #Metrics
                   "Max NonPeak (1)",
                   "Peak:NonPeak (1)",
                   "Max Peak (2)",
                   "Max NonPeak (2)",
                   "Peak:NonPeak (2)",
                   "Shared (PyDegradome)",
                   "Mean Peak Test (Scaled)",
                   "Mean Max Tx Ctrl (Scaled)",
                   "Chr", #Coordinates
                   "Strand",
                   "P. start",
                   "P. stop",
                   "P. Width",
                   "Gene start",
                   "Gene end",
                   "Feat. start",
                   "Feat. end",
                   "CDS len.",
                   "Tx len."
                   )

## Plot settings
sample_color <- c("#FC4E2A", "#E31A1C", "#2c84a8ff", "#214b7fff")

## Highlight color
col_hl <- "#e5ff00bf" #f2ff00cc #e5ff00d9

## Color for bases
base_col <- c(A = "#33A02C", C = "#1F78B4", U = "#E31A1C", G = "#FF7F00")

## font settings
i.font <- ps.options()$family

## minimum width of plot
min_plot_width <- 1000

vars2delete <- c("cds_max",
                 "comp_list_sub",
                 "df_dup",
                 "df_np",
                 "df_np_list",
                 "df_uniq",
                 "dup_ID",
                 "i",
                 "i.cols",
                 "i.id",
                 "i.log.dup",
                 "i.not",
                 "iannot",
                 "iannot_locus_yes",
                 "icoordinates",
                 "idescription",
                 "idir",
                 "idir_prev",
                 "idf",
                 "ilog",
                 "isel",
                 "istr",
                 "mart",
                 "min_tx",
                 ## "mirna",
                 "n_occur_id",
                 "option_list",
                 "opt",
                 "opt_parser",
                 "sel",
                 "temp1",
                 "temp2",
                 "temp_df",
                 "tmp",
                 "tmp2",
                 "treatments",
                 "tx_info",
                 "txdb",
                 "x",
                 "j",
                 "k")
rm(list=vars2delete)
save(list=ls(),
     file=file.path(Rvarious_dir,"Initialization_variables.RData"))
ibase <- env["ibase"]
supp_data_dir <- env["supp_data_dir"]
save(list=c("ibase","supp_data_dir"),
     file=paste0("Env_variables/minimal_variables_", ibase, ".RData"))
