## 01-Annotation-shared_peak_classification.R

## Description: Process PyDegradome results first by matching  coordinates to gene information (Transcript, feature, coordinates). Then, calculates peak signal to background ratio. Finally classifies peaks based on whethere they are found on both replicates (i.e. coordinates overlap). Main output file "Comparison_*"
## The above are performed in 3 for loops. A dictionary-like of reads is built to speed up computation under different settings (prevents reading bigwig files).

## Output dirs:
## - output_02/01-pyDegradome_processed/<sample1_rep1>-<sample2_rep1>/
## - output_02/02-pyDegradome_pooled/
## Output files:
## - <sample1_rep1>-<sample2_rep1>_<conf>_<win>_<mf>_Annotated
## - <sample1_rep1>-<sample2_rep1>_<conf>_<win>_<mf>_Annotated
## - Comparison_<sample1_rep1>-<sample2_rep1>_and_<sample1_rep2>-<sample2_rep2>_<conf>_<win>_<mf>
## Pydegradome settings: conf: Confidence level; win: width of peak search window; mf: Multiplicative factor.

## Examples:
## - output_02/01-pyDegradome_processed/SRR10759112-SRR10759114/SRR10759112-SRR10759114_0_95_4_4_Annotated
## - output_02/02-pyDegradome_pooled/Comparison_SRR10759112-SRR10759114_and_SRR10759113-SRR10759115_0_95_4_4

## Libraries
suppressPackageStartupMessages(library(ChIPpeakAnno)) #For shared peak analysis.  toGRanges function 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gsubfn))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(mgsub))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(PostPyDeg))
library(optparse)

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
##-------------------------------------------------------------------------
##Load files
##Import file with transcript information
txdb <- loadDb(file.path(env["supp_data_dir"],"R/txdb_object"))

annoData.gmodel <- toGRanges(txdb, feature = c("geneModel"))
df.gmodel <- data.frame(annoData.gmodel)

q                                        # Change columns 2,3 names to match function argument
colnames(df.gmodel) <- c("seqnames","feature_start",
                         "feature_end","width","strand",
                         "tx_name","feature_type")

df.gmodel$tx_name <- toupper(df.gmodel$tx_name)



dg_bigwig_all <- list()
for (i.sample in sample_list) {
                                        #Load bigWig [genome]
    bigwig_f <- import(file.path(bigwig_dir,paste0(i.sample,"_G_f.bw")))
    strand(bigwig_f) <- "+"

    bigwig_r <- import(file.path(bigwig_dir,paste0(i.sample,"_G_r.bw")))
    strand(bigwig_r) <- "-"

    dg_bigwig_all[[i.sample]] <- c(bigwig_f,bigwig_r)
}

##Columns to order data frame before exporting
cols.indx <- c("ID",
               "chr",
               "strand",
               "peak_start",
               "peak_stop")

i.cols.tx <- c("feature_start","feature_end","width","tx_name","feature_type")
##--------------------------------------------------
## Add information to peaks, associated transcript and feature
## Output:  'pydeg_annotated_dir'
cat("* Add information to peaks, associated transcript and feature\n")
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    for(i.test in c("test","nc")) {
        for(i.comp in comp_list[[i.test]]) {
                                        # Directory for output
            i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"),c("","-"))
            i.conf_f <- gsub("\\.","_",i.conf)
            out_dir <- file.path(pydeg_processed_dir,i.comp_f)
            
            file.suffix <- paste(i.comp_f,i.conf_f,"4",i.MF,sep="_")
            file.suffix2 <- paste(i.comp,i.conf,"4",i.MF,sep="_")
            ##-----------------------------------
            ##Avoid re-doing annotaion
            out_file <-  file.path(out_dir, paste0(file.suffix,"_Annotated"))
            input_file  <-  file.path(pydeg_dir, file.suffix2)
            ##Test for input and output file
            if(!file.exists(out_file) && file.exists(input_file)) {
                dir.create(out_dir)

                cat("\tProcessing ",i.comp,"...\n")
                cat("\tMF=",i.MF," conf=", i.conf, " Set=", i.test, "\n")
                cat("\t\tLoading", basename(input_file),"\n")

                df  <-  read.table(input_file,
                                   comment.char = "[",
                                   header=TRUE)

                colnames(df)[1] <- "chr"

                                        #Add peak width colmun
                df$peak_width <- df$peak_stop-df$peak_start

                                        #Add annotation
                cat("\t\tAdding annotation...\n")
                ##Find peaks that overlap (within) a feature
                gr.df <- GRanges(#
                    seqnames = df$chr,
                    strand = df$strand,
                    ranges = IRanges(start = df$peak_start,
                                     end = df$peak_stop))
                ##Find peak overlaps (all)
                peak_overlap <- findOverlaps(gr.df,
                                             annoData.gmodel,
                                             type="within",
                                             select="all",
                                             ignore.strand=FALSE)
                
                peak_yes <- !is.na(peak_overlap)

                ##Indices for peaks that overlap with a gene feature
                indx_peak_overlap <- peak_overlap[peak_yes]

                ##Select peaks that mapped to a feature
                df_yes <- df[queryHits(indx_peak_overlap),]
                ## df_no <- df[!peak_yes,]

                ##Select matching gene/feature information
                df.gmodel_sub <- df.gmodel[subjectHits(indx_peak_overlap),]
                df.gmodel_sub <- dplyr::select(df.gmodel_sub ,
                                               all_of(i.cols.tx))
                df.gmodel_sub$ID <- gsub("\\.[0-9]","",
                                         df.gmodel_sub$tx_name)

                ##Combine data frames [colum-wise]
                df_annotated <- cbind(df_yes, df.gmodel_sub)
                df_annotated$indx <- paste0(df_annotated$ID,
                                            df_annotated$peak_start,
                                            df_annotated$peak_stop)
                df_annotated <- df_annotated[with(df_annotated,
                                                  order(indx)),]
                
                ##Import file with transcript information
                tx_info <- fread(file.path(env["supp_data_dir"],
                                           "R/Transcript_information.txt"))

                tx_name_representative <- tx_info[tx_info$rep_gene==1,tx_name]


                ##remove duplicate peaks based on representative gene model
                dfs.ID <- SplitDFbyNpeaks(df_annotated, ref_col = "indx", keep_ref=TRUE)
                df_dup <- dfs.ID$df_dup
                df_uniq <- dfs.ID$df_uniq

                                        # Work on genes with isoforms
                dup_indx <- unique(df_dup$indx)
                df_np_list <- list()
                for(k in seq_along(dup_indx)) {
                    i.indx <- dup_indx[k]
                    temp_df <- df_dup[df_dup$indx==i.indx,]
                    i.log <- temp_df$tx_name %in% tx_name_representative
                    if(sum(i.log)==1) {
                        temp_df <- temp_df[i.log,]
                        temp_df$rep_gene_note <- "Matches"
                    } else {
                                        #pick the first entry and add note
                        temp_df <- temp_df[1,]
                        temp_df$rep_gene_note <- "Doesn't match"
                    }
                    df_np_list[[k]] <- temp_df 
                }
                df_np <- do.call(rbind, df_np_list)

                                        #Work on genes without isoform
                df_uniq$rep_gene_note <- "Matches"

                df_annotated <- rbind(df_uniq,df_np)
                df_annotated <- dplyr::select(df_annotated, -indx)
                
                ##Merge information to df.gmodel
                df_annotated <- merge(#
                    df_annotated,
                    dplyr::select(tx_info,
                                  c("tx_name",
                                    "gene_region_start",
                                    "gene_region_end",
                                    "cds_len",
                                    "tx_width",
                                    "gene_name",
                                    "rep_gene",
                                    "Description")),
                    by="tx_name", all.x =TRUE)

                                        # Remove entries with NAs
                df_annotated <- df_annotated[complete.cases(df_annotated),]
                
                write.table(df_annotated,
                            out_file ,
                            sep="\t",row.names =FALSE)
                cat("\t\tSaved as", basename(out_file),"\n")
            }# Test for input and output
        }#i.comp
    } #i.test
}# Loop over pydeg settings
##---------------------------------------------------
## Create a dataframe with non peak coordinates
## Output: 'pydeg_np_dir'
cat("* Create a dataframe with non peak coordinates\n")
for(i.test in c("test","nc")) {
    for(i.comp in comp_list[[i.test]]) {
        for (i in seq_along(MF_list)) {
            i.MF <- MF_list[i]
            i.conf <- conf_list[i]
            ##--------------------------------------------------
            i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"),c("","-"))
            i.conf_f <- gsub("\\.","_",i.conf)
            file.suffix <- paste(i.comp_f,i.conf_f,"4",i.MF,sep="_")
            
                                        #Read data frame
            out_file <- file.path(pydeg_np_dir,
                                  paste("NP_coordinates", i.comp_f,
                                        i.conf_f,"4", i.MF,sep="_"))
            
            input_file <- file.path(pydeg_processed_dir,i.comp_f,
                                           paste0(file.suffix,"_Annotated"))
            if(!file.exists(out_file) && file.exists(input_file)) {
                ## if(file.exists(input_f)) {
                    cat("\tMF=",i.MF," conf=", i.conf, " Set=", i.test, "\n")
                    df <- read.table(input_file, header=TRUE)
                    
                                        # Select only complete cases and
                                        # those entries with ID                         
                    df.ID <- df[grepl("^AT",df$ID),]

                                        #Sort based on peak_start
                    df.sort <- df.ID[with(df.ID,order(ID,peak_start)),]

                                        #Create a indx vector
                    df.sort$indx_dup <-  paste0(df.sort$ID,
                                                df.sort$chr,
                                                df.sort$strand,
                                                df.sort$gene_region_start,
                                                df.sort$gene_region_end)

                                        #Count and select duplicated entries
                    dfs.ID <- SplitDFbyNpeaks(df.sort,"indx_dup")
                    df_dup <- dfs.ID$df_dup
                    ## df_uniq <- dfs.ID$df_uniq
                    
                    dup_ID <- unique(df_dup$tx_name)
                    df_np_list <- list()
                    if (length(dup_ID)>0) {
                        for(k in seq_along(dup_ID)) {
                            i.id <- dup_ID[k]
                            temp_df <- df_dup[df_dup$tx_name==i.id,]
                            g.start <- temp_df$gene_region_start[1]
                            g.end <- temp_df$gene_region_end[1]
                            i.start <- temp_df$peak_stop+1
                            i.end <- temp_df$peak_start-1
                            if(i.end[1]<g.start) {
                                g.start <- NULL
                                i.end <- i.end[-1]
                            }
                            if(i.start[length(i.start)]>g.end) {
                                g.end <- NULL
                                i.start <- i.start[-length(i.start)]
                            }
                            np.start <- c(g.start,i.start)
                            np.stop <- c(i.end,g.end)
                            i.chr <- temp_df$chr[1]
                            i.strand <- temp_df$strand[1]
                            i.ID <- temp_df$tx_name[1]
                            tmp_df_np_coor <- data.table(chr=i.chr,
                                                         strand=i.strand,
                                                         ID=i.id,
                                                         non_peak_start=np.start,
                                                         non_peak_stop=np.stop)

                                        #Remove negative ranges
                            i.diff <- tmp_df_np_coor$non_peak_stop - tmp_df_np_coor$non_peak_start
                            i.log <- i.diff<0
                            tmp_df_np_coor <- tmp_df_np_coor[!i.log,]

                            df_np_list[[k]] <- tmp_df_np_coor
                        } #Loop over ID
                        df_np <- do.call(rbind, df_np_list)
                        write.table(df_np, out_file,
                                    sep="\t", row.names=FALSE, quote=FALSE )
                        cat("\t\tSaved as", basename(out_file),"\n")
                    }#if len dup_ID > 0
            } # Check for input and output
        } # Loop over pydeg settings
    }#i.comp (test and ctl)
}#i.test
##--------------------------------------------------
## Calculate ratio between peak reads and reads outside the peak region (signal:noise)
## Creates a dictionary-like of coordinates and reads to speed up computation
cat("* Calculate peak to non-peak ratio \n")
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    for(i.test in c("test","nc")) {
        for(i.comp in comp_list[[i.test]]) {
            i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"), c("","-"))
            i.conf_f <- gsub("\\.","_",i.conf)
            out_dir <- file.path(pydeg_processed_dir,i.comp_f)
            file.suffix <- paste(i.comp_f,i.conf_f,"4",i.MF,sep="_")
            out_file <- file.path(out_dir, paste0(file.suffix,"_Processed"))
            input_file <- file.path(out_dir,paste0(file.suffix,"_Annotated"))

            ## Test if pydeg_annotated exist
            if (!file.exists(out_file) && file.exists(input_file)) {
                                    
                    cat("\tMF=",i.MF," conf=", i.conf, " Set=", i.test, "\n")
                    cat("\t\tLoading", basename(input_file),"\n")
                    cat("\t\tProcessing ",i.comp,"...\n")
                    df <- read.table(input_file,
                                     sep = "\t",
                                     comment.char = "",
                                     row.names = NULL,
                                     header = TRUE)
                    
                    cat("\t\tCheck Ref file & Add max read in peak\n")
                    df <- checkPeakRef(df,i.comp)

                    ## Add S/N ratio
                    cat("\t\tAdd S/N ratio \n")
                    ##--------------------------------------------------
                    cat("Check ref & Add highest signal outside of the peak region\n")
                    ##--------------------------------------------------
                                        #Subset entries with ID
                    df.ID <- df[grepl("^AT",df$ID),]

                                        #Create indx_dup to select duplicates
                    df.ID$indx_dup <-  paste0(df.ID$ID,
                                              df.ID$chr,
                                              df.ID$strand,
                                              df.ID$gene_region_start,
                                              df.ID$gene_region_end)

                    dfs.ID <- SplitDFbyNpeaks(df.ID,"indx_dup")
                    df_dup <- dfs.ID$df_dup
                    df_uniq <- dfs.ID$df_uniq
                    ## if(nrow(df_dup) > 0)
                    if(!is.null(df_dup) && nrow(df_dup) > 0) {
                        ##==================================================
                        cat("Calculating max signal outside peak  (Peaks 2+)...\n")

                                        # Import bigwig data based on i.comp variable
                        sample_name <- SampleName(i.comp)
                        bigwig <- dg_bigwig_all[[sample_name]]
                        ##--------------------------------------------------
                                        #Calculate read outside peak
                                        #Region-wise
                        cat("\tCalculating max signals outside peak (Region-wise)...\n")
                                        #Load coordinate file
                        pydeg_np_region_f <- file.path(pydeg_np_dir,
                                                       paste("NP_coordinates",
                                                             file.suffix,sep="_"))
                        cat("\t\tLoading",
                            basename(pydeg_np_region_f),"\n")
                        pydeg_np_region <- read.table(pydeg_np_region_f,
                                                      sep = "\t",
                                                      row.names = NULL,
                                                      header = TRUE)
                                        #Add non peak region
                        pydeg_np_region <- addMaxNonPeakRegion(#
                            df=pydeg_np_region,
                            bigwig=bigwig,
                            i.comp=i.comp,
                            core=env$core)
                        cat("\tFinished calculating max signals outside peak (Region-wise)\n")
                        ##--------------------------------------------------
                                        #Calculate read outside peak
                                        #Gene-wise
                        cat("\tCalculating max signals outside peak (Gene-wise)...\n")
                                        #Vector of unique IDs to loop
                        dup_ID <- unique(df_dup$tx_name) #ID
                        pydeg_np_gene <- addMaxNonPeakGene(pydeg_np_region, dup_ID, core=env$core)
                        pydeg_np_gene <- pydeg_np_gene[!is.na(pydeg_np_gene$ID),]
                        ##--------------------------------------------------
                                        #Remove i.cols and bigwig
                        rm(bigwig)
                        cat("\tFinished calculating max signal outside peak (Peaks 2+)\n")
                        ##==================================================
                        cat("\tAdding max_non_peak_gene to original data frame (df_dup)...\n")
                                        #Remove duplicates
                        pydeg_np_gene <- unique(pydeg_np_gene)
                        df_dup <- MergeDFs(dfA=df_dup,
                                           dfB=pydeg_np_gene,
                                           ref_col="max_np_gene")
                        ## df_dup <- dplyr::select(df_dup,-indx_dup)
                        cat("Finish adding max_non_peak_gene to original data frame\n")
                    } else {
                        df_dup <- NULL
                    }
                    ##----------------------------------------
                    cat("\t\tCalculating Max Non Peak (Peaks 1)...\n")
                    if(nrow(df_uniq) > 0) {
                        df_uniq <- checkNonPeakRefOne(df_uniq,i.comp,cols.indx)
                        ## df_uniq <- dplyr::select(df_uniq,-indx_dup)
                    }else {
                        df_uniq <- NULL
                    }
                    cat("\t\tFinished calculating Max Non Peak (Peaks 1)\n")
                    ##----------------------------------------

                    df.NonID <- df[!grepl("^AT",df$ID),]#!is.na(df$ID)
                    if(nrow(df.NonID) > 0) {
                        df.NonID$max_np_gene <- NA
                    }else {
                        df.NonID <- NULL
                    }
                    cat("\t\tMerging dataframes...\n")

                    df <- rbind(df_dup,df_uniq,df.NonID)

                                        #Change np_gene if 0 to 1 to avoid inf
                    df <- as.data.table(df)
                    df[max_np_gene == 0,max_np_gene := 1 ]
                    ##--------------------------------------------------
                    ##Calculate max_non_peak_ratio
                    cat("Calculate max_non_peak_ratio\n")
                    df$max_non_peak_ratio <- df$max_peak / df$max_np_gene

                    ## Save table
                    cat("Save table \n")
                    write.table(df,out_file,
                                sep="\t",row.names = FALSE )
                    cat("\t\tSaved as", basename(out_file),"\n")
            }# Test for input and output file
        } # i.comp
    } # i.test
} # Loop over pydeg settings
##--------------------------------------------------
## Generate shared peak dataset
## For peaks found in only one sample, it calculates the reads (on the missing sample) using the same coordinates .
## These are latter treated as a peak, though technically no significant difference was found by pyDegradome
cat("* Generate shared peak dataset\n")
peak_num_df <- data.frame() #Stacked tables of peak width
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    pydeg_list <- list() ##List of pyDegradome output
    peak_range_list <- list() ##List of pyDegradome peaks (GRanges)
    i.conf_f <- gsub("\\.","_",i.conf)
    for(i.test in c("test","nc")) {
        cat("\tMF=",i.MF," conf=", i.conf, " Set=", i.test, "\n")
        ##Generate lists from pyDegradome output for comparing pairs
        cat("\t**Generate lists from pyDegradome output for comparing pairs")
        for(i.comp in comp_list[[i.test]]) {
            i.comp_f <- mgsub::mgsub(i.comp,c("t_","_c_"),c("","-"))
            i.conf_f <- gsub("\\.","_",i.conf)
            file.suffix <- paste(i.comp_f,i.conf_f,"4",i.MF,sep="_")
            out_dir <- file.path(pydeg_processed_dir,i.comp_f)
            input_file <- file.path(out_dir,
                                           paste0(file.suffix,"_Processed"))

                                        #Test if pydeg_processed_f exists
            if(file.exists(input_file)) {
                df <- fread(input_file)
                
                ##Peak number
                cat("Peak number \n")
                peak_num_df <- rbind(peak_num_df,
                                     data.frame(count=nrow(df),
                                                sig=i.conf,
                                                mf=i.MF,
                                                test_type=i.test,
                                                comp_pair=i.comp))
                
                ##Convert the results into GRanges obj
                cat("Convert the results into GRanges obj \n")
                peak_range_list[[i.comp]] <- GRanges(#
                    seqnames = df$tx_name,
                    strand = df$strand,
                    ranges = IRanges(#
                        start=df$peak_start,
                        end = df$peak_stop))
                ##Keep raw input + width, annotation) as data.frame for filtering later.
                cat("Keep raw input + width \n")
                pydeg_list[[i.comp]] <- df
            }# if pydeg_processed_f exists
        } #i.comp
        ##[end] Calculate peak width
        ##--------------------------------------------------
        names(peak_range_list) <- mgsub::mgsub(names(peak_range_list),c("t_","_c_"),c("","-"))
        ##--------------------------------------------------

        ##Filter pyDegradome results by shared peak between replicates
        cat("\t**Filter pyDegradome results by shared peak \n")
        ## Get unique combinations of pairs
        tmp1 <- mgsub::mgsub(comp_list[[i.test]], c("t_","_c_"), c("","-"))
        pairs_comb <- expand.grid(tmp1[1], tmp1[2:length(tmp1)])
        tmp2 <- paste(pairs_comb[,1], pairs_comb[,2], sep="-")
        isamples <- strsplit(tmp2, "-")
        uniq_indx <- sapply(isamples, function(i) length(unique(i)) == 4)
        pairs_comb_uniq <- pairs_comb[uniq_indx,]
        
        for(i.pair in seq_len(nrow(pairs_comb_uniq))) {
            i_pair1 <- as.character(pairs_comb_uniq[i.pair,1])
            i_pair2 <- as.character(pairs_comb_uniq[i.pair,2])
            out_file <- file.path(pydeg_pooled_dir,
                                  paste("Comparison",i_pair1,"and",
                                        i_pair2,i.conf_f,"4",i.MF,sep="_"))
            
            if(!file.exists(out_file)) {

                cat("\t\tGenerating, ",
                    basename(out_file),  "\n")
                i.file1 <- paste(file.path(pydeg_processed_dir, i_pair1, i_pair1),
                                 i.conf_f,"4",i.MF,"Processed",sep = "_")
                
                i.file2 <- paste(file.path(pydeg_processed_dir, i_pair2, i_pair2),
                                 i.conf_f,"4",i.MF,"Processed",sep = "_")
                
                if(file.exists(i.file1) && file.exists(i.file2)) {

                    cat("\t\tFinding overlaps \n")
                                        #Get overlapped peaks between 2 samples
                    overlaps <- findOverlaps(#
                        peak_range_list[[i_pair1]],
                        peak_range_list[[i_pair2]],
                        ignore.strand=FALSE)
                    
                    if (!is.null(overlaps) && length(overlaps)>0) {
                        icols <- c("max_np_gene","max_peak", "max_non_peak_ratio","max_count_test")
                                        #Comparisons tested
                        i.comp_list<- list(#
                            i.comp1= i_pair1,
                            i.comp2= i_pair2
                        )
                        
                                        #Indices of overlap peaks
                        sel_list <- list(
                            sel1=queryHits(overlaps),
                            sel2= subjectHits(overlaps)
                        )
                        ##--------------------------------------------------
                        pydegShared_peak_list <- list()
                        pydegSingleRep_peak_list <- list()
                        for (j in seq_along(sel_list)) {
                            ## Define variables for indexing
                            j.not <- seq_along(sel_list)[!seq_along(sel_list) %in% j]
                            i.comp <- gsub("-","_c_",i.comp_list[[j]]) %>% paste0("t_",.) 
                            i.comp.other <- gsub("-","_c_",i.comp_list[[j.not]]) %>% paste0("t_",.)

                            ## Indices of overlaps
                            sel <- sel_list[[j]]

                            ##Change column names
                            pydeg <- data.frame(pydeg_list[[i.comp]])
                            setnames(pydeg,
                                     old=icols,
                                     new=paste(icols,j,sep="_"))
                            
                            
                            ## Work on peaks found on both samples
                            pydegShared_peak <- pydeg[sel,]
                            pydegShared_peak_list[[j]] <- pydegShared_peak
                            
                            ## Work on peaks found in one sample
                            pydegSingleRep_peak <- pydeg[-sel,]
                            pydegSingleRep_peak$shared <- j
                            
                            rm(pydeg)
                            ##select entries with ID
                            df.ID <- pydegSingleRep_peak[grepl("^AT",
                                                               pydegSingleRep_peak$ID),]
                            
                            ##Create indx_dup to select duplicates
                            df.ID <- within(df.ID,
                                            indx_dup <- paste(ID,chr,strand,gene_region_start,gene_region_end,
                                                              sep=""))
                            dfs.ID <- SplitDFbyNpeaks(df.ID,"indx_dup")
                            df_dup <- dfs.ID$df_dup
                            df_uniq <- dfs.ID$df_uniq
                            
                            ##Check that the dataframe is not empty
                            if(!is.null(df_dup) && nrow(df_dup) > 0) {

                                cat("\t\tCalculating max signal outside peak  (Peaks 2+)...\n")

                                        # Import bigwig data based on i.comp variable
                                sample_name <- SampleName(i.comp.other)
                                bigwig <- dg_bigwig_all[[sample_name]]
                                        #Calculate read outside peak
                                        #Region-wise
                                cat("\t\tCalculating max signals outside peak (Region-wise)...\n")
                                        #Set i.comp to sample 2 to load coordinates of non-peaks
                                pydeg_np_region_f <- file.path(#
                                    pydeg_np_dir,
                                    paste("NP_coordinates",
                                          file.suffix,sep="_"))
                                cat("\t\tLoading",
                                    basename(pydeg_np_region_f),"\n")
                                
                                pydeg_np_region <- read.table(#
                                    pydeg_np_region_f,
                                    sep = "\t",
                                    row.names = NULL,
                                    header = TRUE)

                                        #Set i.comp to sample 1 to load reference file
                                pydeg_np_region <- addMaxNonPeakRegion(df=pydeg_np_region,
                                                                       bigwig=bigwig,
                                                                       i.comp=i.comp.other,
                                                                       core=env$core)
                                cat("\t\tFinished calculating max signals outside peak (Region-wise)\n")
                                        #Vector of unique IDs to loop
                                dup_ID <- unique(df_dup$tx_name) #ID
                                        #Calculate read outside peak
                                        #Gene-wise
                                cat("\t\tCalculating max signals outside peak (Gene-wise)...\n")
                                pydeg_np_gene <- addMaxNonPeakGene(#
                                    pydeg_np_region, dup_ID, core=env$core)
                                ##--------------------------------------------------
                                        #Remove i.cols and bigwig
                                rm(bigwig)
                                cat("\t\tFinished calculating max signal outside peak  (Peaks 2+)\n")
                                ##==================================================
                                cat("\t\tAdding max_non_peak_gene to original data frame (df_dup)...\n")
                                        #Remove duplicates
                                df_np_gene <- unique(pydeg_np_gene)
                                df_dup <- MergeDFs(dfA=df_dup,
                                                   dfB=df_np_gene,
                                                   ref_col="max_np_gene")
                                cat("\t\tFinished adding max_non_peak_gene to original data frame\n")
                            } else {
                                df_dup <- NULL
                            }
                            ##Working on genes with a single peak
                            cat("\t\tCalculating Max Non Peak (Peaks 1)...\n")
                            if(!is.null(df_uniq) && nrow(df_uniq) > 0) {
                                df_uniq <- checkNonPeakRefOne(#
                                    df_uniq,
                                    i.comp=i.comp.other,
                                    cols.indx)
                                ## df_uniq <- dplyr::select(df_uniq,-indx_dup)
                            } else {
                                df_uniq <- NULL
                            }
                            cat("\t\tFinished calculating Max Non Peak (Peaks 1)\n")
                            ##==================================================
                            ##Merging all the dataframes in to the original df
                            cat("\t\tMerging dataframes...\n")
                            pydegSingleRep_peak  <- rbind(df_dup,df_uniq)
                            
                            ##Change np_gene if 0 to 1 to avoid inf
                            setDT(pydegSingleRep_peak)
                            pydegSingleRep_peak[max_np_gene == 0,
                                                max_np_gene := 1 ]
                            setDF(pydegSingleRep_peak)
                            ##==================================================
                            
                            ##Calculate max peak for the opposite replicate
                            cat("\t\tCheck Ref & Calculate max peak for rep",
                                j.not, " (", i.comp.other ,")\n")
                            
                            pydegSingleRep_peak <- checkPeakRef(df=pydegSingleRep_peak, i.comp=i.comp.other)
                            
                            ##Calculate max_non_peak_ratio
                            pydegSingleRep_peak$max_non_peak_ratio <- pydegSingleRep_peak$max_peak / pydegSingleRep_peak$max_np_gene

                            ##Fill max_count_test with NA
                            pydegSingleRep_peak$max_count_test <- NA

                            ##rename columns
                            setnames(pydegSingleRep_peak,
                                     old=icols,
                                     new=paste(icols,j.not,sep="_"))
                            
                            ##Add df to list
                            pydegSingleRep_peak_list[[i]] <- pydegSingleRep_peak
                        }
                        pydeg_single_peak <- do.call(rbind,
                                                     pydegSingleRep_peak_list)
                        
                        ##--------------------------------------------------
                        ## Work on peaks found in both samples
                        pydeg_shared_peak <- cbind(#
                            pydegShared_peak_list[[1]],
                            dplyr::select(pydegShared_peak_list[[2]],
                                          all_of(paste(icols,"2",sep="_"))))
                        pydeg_shared_peak$shared <- 3

                                        # width takes average of both samples
                        pydeg_shared_peak$width <- rowMeans(cbind(#
                            pydegShared_peak_list[[1]][["width"]],
                            pydegShared_peak_list[[2]][["width"]]))

                                        # remove indx column
                        isel <- colnames(pydeg_shared_peak) %in%
                            colnames(pydeg_single_peak)
                        pydeg_shared_peak <- pydeg_shared_peak[,isel]

                        cat("\t\tCombining data frames... \n")
                        pydeg_all_peak<- rbind(pydeg_single_peak,
                                               pydeg_shared_peak)

                        ##Save processed pyDegradome result
                        cat("Save processed pyDegradome result \n")
                        pair_comparison <- paste(pairs_comb_uniq[i.pair,1],
                                                 "and",
                                                 i_pair2,sep = "_")
                        write.table(pydeg_all_peak,
                                    out_file, sep="\t",row.names =FALSE)
                        ##Peak number
                        cat("\t\tPeak number \n")
                        peak_num_df <- rbind(peak_num_df,
                                             data.frame(count=nrow(pydeg_shared_peak),
                                                        sig=i.conf,mf=i.MF,test_type=i.test,
                                                        comp_pair=pair_comparison))
                    }#If there are overlaps
                }# if both files exists pair_df 1&2
            }#if pydeg_pairs_comparison_f exists   
        }#i.pair
    } #i.test
} # Loop over pydegradome settings
##--------------------------------------------------
##Save peak number table.
peak_num_df <- unique(peak_num_df)
write.table(peak_num_df,file.path(summary_dir,"Annotated_peak_number"),
            sep="\t", row.names=FALSE, quote=FALSE )
## }#Loop over comp_list__sub
