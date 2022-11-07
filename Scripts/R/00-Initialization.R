## Libraries
suppressPackageStartupMessages(library(here))
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
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Working directory", metavar="character"),
    make_option(c("-b", "--base"), type="character", default=NULL, 
                help="Project basename", metavar="character")
) 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$wd)
## Path settings

env.raw <- read.table(here(paste0("Env_variables/Degradome_",opt$base,".txt")),
                      sep = "=", quote="")
env.raw[,2] <- gsub('" "', '"\\, "', env.raw[,2])
isel <- grepl("\\(",env.raw[,2])
env.raw1 <- env.raw[isel,]
env.raw2 <- env.raw[!isel,]

for (i in seq_len(nrow(env.raw1))) {
    istr <- env.raw1[i,2]
    assign(env.raw1[i,1],
           eval(parse(text=paste0("c",gsub('" "', '"\\, "', istr)))))
  }

env <- str_trim(env.raw2$V2)
names(env) <- str_trim(env.raw2$V1)
env <- as.list(env)



## If `core` is not defined in env, use a single core.
if(!is.element("cores",names(env))) {
    env$cores <- 1
} else {
    env$cores <- as.numeric(env$cores)
}

## Input path
#Degradome BigWig
bam_dir <- file.path(env["output_dirB"], "04-bam_genomic_tophat")
bam_tr_dir <- file.path(env["output_dirB"], "06-bam_transcript_tophat")
bigwig_dir <- file.path(env["output_dirB"],"04-bigwig_genomic_tophat")
pydeg_dir <- file.path(env["output_dirB"],"05-pyDegradome_tophat")
htseq_dir <- file.path(env["output_dirB"], "07-htseq_genomic_tophat")


## Output path
pydeg_processed_dir <- file.path(env["output_dirR"],"01-pyDegradome_processed")
pydeg_pooled_dir <- file.path(env["output_dirR"],"02-pyDegradome_pooled")
report_dir <- file.path(env["output_dirR"],"03-Report")
summary_dir <- file.path(report_dir,"Summary")

## Peak reads
pydeg_reads_dir <- file.path(env["supp_data_dir"],"Degradome_reads")

## Coordinates of regions with no peaks
pydeg_np_dir <- file.path(env["supp_data_dir"],"pyDegradome_NonPeaks_coordinates")

## Max peak reads 
maxP_dir <- file.path(env["supp_data_dir"], "pyDegradome_maxPeak_CPM")

## Max transcript reads
maxR_dir <- file.path(env["supp_data_dir"], "pyDegradome_maxTx_CPM")

## Root dir to store degradome plots
pydeg_dplot_dir <- file.path(report_dir, "Dplots")

## Fasta dir
peakSeq_dir <- file.path(env["supp_data_dir"], "Peak_sequences")

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

## Settings used for pyDegradome analysis
## MF_list <- 3:4 # Multiplicable factor
## conf_list <- c(0.95,0.99) # Confidence level
pydeg_settings <- unlist(strsplit(pydeg_script_settings," "))
MF_list <- as.numeric(pydeg_settings[rep(c(FALSE,TRUE),
                                         length(pydeg_script_settings))])
conf_list <- as.numeric(pydeg_settings[rep(c(TRUE,FALSE),
                                         length(pydeg_script_settings))])
## Sample name list
## control_samples <- c("SRR10759114","SRR10759115")
names(control_samples) <- control_samples_name
## test_samples <- c("SRR10759112","SRR10759113")
names(test_samples) <- test_samples_name

idf <- expand.grid(test=test_samples,nc=control_samples)
comp_list <- list()
for (i in seq_along(colnames(idf))) {
  i.not <- seq_along(colnames(idf))[!seq_along(colnames(idf)) %in% i]
  i.comp <- colnames(idf)[i]
  i.comp_other <- colnames(idf)[i.not]
  comp_list_sub <- NULL
  for (j in seq_len(nrow(idf))) {
    x <- paste("t",idf[j,i.comp],"c",idf[j,i.comp_other],sep="_")
    comp_list_sub <- c(comp_list_sub,x)
  }
  comp_list[[i.comp]] <- comp_list_sub
}

tmp <- data.frame(sapply(comp_list,c))
tmp2 <- list(as.data.frame(combn(tmp[["test"]],2)),
             as.data.frame(combn(tmp[["nc"]],2)))
comp_pair_list <- list()
for (k in seq_along(tmp2[[1]])) {
    comp_pair_list[[k]] <- list(test=paste(tmp2[[1]][,k],collapse ="_and_"),
                                nc=paste(tmp2[[2]][,k],collapse ="_and_"))
}

treatments <- gsub(" \\[[0-9]\\]","",
                   c(names(control_samples),names(test_samples)))

idesign <- data.frame(sample=c(control_samples,test_samples),
                      replicate=c(seq_along(control_samples),seq_along(test_samples)),
                     treatment=treatments)
rownames(idesign) <- seq_len(nrow(idesign))

sample_list <- c(test_samples,control_samples)

dg_bigwig_list <- structure(file.path(bam_dir,sample_list), names = sample_list)

dg_bam_list <- structure(file.path(bam_dir,paste0(sample_list,".bam")),
                         names = sample_list)

dg_bam_tr_list <- structure(file.path(bam_tr_dir,paste0(sample_list,"_uni_sort.bam")),
                         names = sample_list)
##--------------------------------------------------
##Build a transcript info data frame
##Data frame with feature, start, stop, gene name information from biomaRt
tx_info_f <- file.path(Rvarious_dir,"Transcript_information.txt")
if (!file.exists(tx_info_f)) {
  mart <- useMart(dataset="athaliana_eg_gene",
               biomart="plants_mart",
               host="plants.ensembl.org")

  idescription <- getBM(
      mart = mart,
      attributes = c(
          'ensembl_transcript_id',
          'transcript_start',
          'transcript_end',
          'external_gene_name',
          'entrezgene_description',
          'interpro_short_description')
  )

  icoordinates <- getBM(
      mart = mart,
      attributes = c('ensembl_transcript_id',
                     'cds_length'))

  iannot <- merge(idescription,icoordinates,
                  by='ensembl_transcript_id',all.x=TRUE)

  colnames(iannot) <- c("tx_name", "gene_region_start", "gene_region_end",
                        "gene_name", "entrezgene_description",
                        "interpro_short_description", "cds_len")
  
  i.cols <- c("ID", "tx_name", "gene_region_start",
              "gene_region_end", "cds_len", "gene_name", "Description")

  iannot_locus_yes <- iannot[iannot$tx_name!="",]
  iannot_locus_yes$ID <- gsub("\\.[0-9]","",iannot_locus_yes$tx_name)

                                        # Create a Description column based on the entrezgene and interpro descriptions
  ilog <- iannot_locus_yes$entrezgene_description==""
  temp1 <- iannot_locus_yes[ilog,]
  temp2 <- iannot_locus_yes[!ilog,]

  temp1$Description <- temp1$interpro_short_description
  temp2$Description <- temp2$entrezgene_description

  temp1 <- dplyr::select(temp1,all_of(i.cols))
  temp2 <- dplyr::select(temp2,all_of(i.cols))
  tx_info <- rbind(temp1,temp2)
  
  
                                        #Create an index column to remove duplicated entries
  tx_info$indx <-  paste0(tx_info$tx_name,
                          tx_info$gene_region_start,
                          tx_info$gene_region_end)
  
  sel <- duplicated(tx_info$indx)
  tx_info <- dplyr::select(tx_info[!sel,],all_of(i.cols))


  ## tx_info <- tx_info[with(tx_info,order(tx_name)),]

  tx_info$tx_width <- tx_info$gene_region_end - tx_info$gene_region_start

                                        #Identify genes with and without isoforms
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
      temp_df <- df_dup[df_dup$ID==i.id,]
      ## tx_max <- max(temp_df$tx_width)
      cds_max <- max(temp_df$cds_len)
      temp_df$rep_gene <- if_else(temp_df$cds_len == cds_max, 1, 0)

      if(length(which(temp_df$rep_gene==1))>1) {
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
  tx_info <- rbind(df_np,df_uniq)
  tx_info <- tx_info[with(tx_info,order(tx_name)),]
  
  fwrite(tx_info,tx_info_f, quote = FALSE, sep = "\t")
}
##--------------------------------------------------
#Make a TxDb object from transcript annotations available as a GFF3 or GTF file
txdb_f <- file.path(Rvarious_dir,"txdb_object")
if (!file.exists(txdb_f)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(#
                                 file = env[["ref_gtf"]], format = "gtf")
    saveDb(txdb,txdb_f)
}

## Import and process miRNA target database
mirna_f <- file.path(Rvarious_dir,"miRNA")
if (!file.exists(mirna_f)) {
    mirna <- read.csv(env$mirna_list)
                                        #Get first record with priority by the evidences.
    mirna$target_type <- factor(mirna$target_type, levels = c("VALIDATED", "PREDICTED", "predicted"))
    mirna <- mirna[order(mirna$target_type),]
    mirna <- mirna[!duplicated(mirna$target_acc),]
    data.table::setnames(mirna,"target_acc","ID")
    saveRDS(mirna,mirna_f)
}

## Generate a sqlite database
ensemblDB_file <- file.path(Rvarious_dir,
                            gsub("\\.gtf",".sqlite",
                                              basename(env[["ref_gtf"]])))
if (!file.exists(ensemblDB_file)) {
    DB <- ensDbFromGtf(gtf=env[["ref_gtf"]], path=Rvarious_dir)
} else {
    DB <- ensemblDB_file
}


## Variables
## At_genome_file <- basename(env[["At_genome"]])
## At_transcript_path <- env[["At_transcript"]]
## ref_gtf <- env[["ref_gtf"]]

#Columns to round up
cols2round <- c("max_non_peak_ratio", "max_non_peak_ratio_2", "max.peak.cpm", "max.read.tx.cpm", "ratioPTx")

## Link to TAIR
tair.prefix <- "<a  target=_blank href=https://www.arabidopsis.org/servlets/TairObject?type=locus&name="
tair.suffix <- ">TAIR</a>"

## Link to PARE
pare.prefix <- "<a  target=_blank href=https://mpss.danforthcenter.org/web/php/pages/GeneAnalysis.php?SITE=at_pare&featureName="
pare.suffix <- "&model=1>"

## Link to Degradome plots
dplot.prefix <- "<a  target=_blank href="

## List of columns for rounding values
cols2round <- c("Peak:NonPeak (1)", "Peak:NonPeak (2)", "Mean Peak Test (CPM)", "Mean Max Tx Ctrl (CPM)", "Peak(T):Max(C)")

## Columns for exporting table of candidates
i.cols.sort <- c("tx_name",
                 "Gene_plot",
                 "Peak_plot",
                 "rep_gene",
                 "feature_type",
                 "gene_name",
                 "ratioPTx",
                 "category_1",
                 "category_2",
                 "MorePeaks",
                 "miRNA",
                 "Description",
                 "comparison",
                 "max_peak_1",#Metrics
                 "max_np_gene_1",
                 "max_non_peak_ratio_1",
                 "max_peak_2",
                 "max_np_gene_2",
                 "max_non_peak_ratio_2",
                 "shared",
                 "max_peak_cpm",
                 "max_read_tx_cpm",
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

i.cols.export <- c("Transcript",#Main
                   "Gene\nplot",
                   "Peak\nplot",
                   "Is rep.",
                   "Feature",
                   "Gene name",
                   "Peak(T):Max(C)",
                   "Cat. 1",
                   "Cat. 2",
                   "Peaks 2>",
                   "miRNA",
                   "Description",
                   "Comparison",
                   "Max Peak (1)", #Metrics
                   "Max NonPeak (1)",
                   "Peak:NonPeak (1)",
                   "Max Peak (2)",
                   "Max NonPeak (2)",
                   "Peak:NonPeak (2)", 
                   "Shared (pyDegradome)",
                   "Mean Peak Test (CPM)",
                   "Mean Max Tx Ctrl (CPM)",
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
col_hl <- "#e5ff00bf" #f2ff00cc#e5ff00d9

## Color for bases
base_col <- c( A="#33A02C",C="#1F78B4",U="#E31A1C", G="#FF7F00")

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
                 "idf",
                 "ilog",
                 "isel",
                 "istr",
                 "mart",
                 "min_tx",
                 "mirna",
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
