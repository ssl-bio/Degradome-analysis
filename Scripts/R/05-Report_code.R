## 05-Report_code.R

## Description: Code to format the output from the Degradome analysis. Code chunks that will be processed by the 05-Report.Rmd file are delimited by ## ---- <name> ---- lines where <name> is the name of the chuck called within the Rmarkdown file.

## ---- setup ----
knitr::opts_chunk$set(echo = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      figure.width=14,
                      results='asis',
                      fig.pos = 'center')

## Set variables
i.conf_f <- gsub("\\.","_",params$i.conf)
i.MF <- params$i.MF

# Set variables
icon_size <- "25px"
pydeg_input_f <- file.path(pydeg_pooled_dir,
                           paste("Pooled",
                                 i.conf_f, "4", i.MF, sep="_"))
##Columns for subset main df (tabs)
cols.gen <- c(i.cols.export[1:7],i.cols.export[22],"PARE",
              i.cols.export[8:9],"TAIR",i.cols.export[10:11])
cols.metrics <- c(i.cols.export[1:4],
                  i.cols.export[c(12:20,11)])
cols.coordinates <- c(i.cols.export[1:4],
                      i.cols.export[c(23:31,11)])
								 

## Import data
pydeg_all <- fread(pydeg_input_f)
setDT(pydeg_all)

## Add Index column for miRNA matching
pydeg_all$Index <- paste0(pydeg_all$comparison,
                         gsub("\\.", "_", pydeg_all$tx_name),
                         pydeg_all$peak_start,"-",
                         pydeg_all$peak_stop)

## Set comparison as factor to loop
pydeg_all$comparison <- as.factor(pydeg_all$comparison)
pydeg_dt_list <- list()
for (i.comparison in levels(pydeg_all$comparison)) {
    pydeg_sub <- pydeg_all[comparison==i.comparison,]

   ## Substitute comparison names
  if(!is.null(pydeg_sub) && nrow(pydeg_sub)>0) {
    i.samples <- unlist(strsplit(gsub("_and_","-",i.comparison),
                                 split="-"))
    i.samples2 <- NULL
    for (i in i.samples) {
      i.sel <-  sample_list %in% i
      names(i) <- names(sample_list)[i.sel]
      i.samples2 <- c(i.samples2,i)
    }
    i.comp <- mergeVector(i.samples2, single=TRUE)
    
    ## Create a data frame of all combination of categories
    i.cat1 <- factor(pydeg_sub$category_1,levels=c(1:4,0))
    i.cat2 <- as.factor(pydeg_sub$category_2)
    cat_df <- expand.grid(levels(i.cat1),
                          levels(i.cat2))
    
    for (j.row in seq_len(nrow(cat_df))) {
        pydeg_cat <- pydeg_sub[pydeg_sub$category_1 == cat_df[j.row,"Var1"],]
        pydeg_cat <- pydeg_cat[pydeg_cat$category_2 == cat_df[j.row,"Var2"],]

        ## Check that the data frame is not empty
        if (!is.null(pydeg_cat) && nrow(pydeg_cat)>0) {
            cat_label <- paste(cat_df[j.row,"Var1"],
                               cat_df[j.row,"Var2"],sep = "-")
            
            ## Limit the number of peaks to show to 5 unless
            ## it is the main combination Cat1=1 and Cat2=A
            if(unique(pydeg_cat$category_1)=="1" &&
               unique(pydeg_cat$category_2)=="A") {
                Nplots <- length(pydeg_cat$category_1)
            } else {
                Nplots <- 5
            }

            ##Sort dataframe by txRatio
            pydeg_cat<- pydeg_cat[with(pydeg_cat,
                                          order(category_1,
                                                category_2,
                                                -ratioPTx)),]
            
            if (nrow(pydeg_cat) < Nplots) {
                np <- nrow(pydeg_cat)
            } else {
                np <- Nplots
            }

            pydeg_cat <- pydeg_cat[1:np,]
            
            ## Preffix and suffix for the plot path
            plot_dir_prefix <- paste(sprintf("%02d",
                                                  as.numeric(row.names(pydeg_cat))),
                                     sep = "_")
            plot_dir_suffix <- paste0(paste(cat_label,
                                            gsub("\\.","_",
                                                 as.vector(pydeg_cat[, "tx_name"])[[1]]),
                                            i.conf_f, "4", i.MF, sep = "_"),
                                      ".pdf")

            ## Plot path
            gene_plot_location <- 
                paste0("Dplots/Gene_",i.comparison,"/",plot_dir_prefix,"_Gene_",plot_dir_suffix)

            peak_plot_location <-
                paste0("Dplots/Peak_",i.comparison,"/",plot_dir_prefix,"_Peak_",plot_dir_suffix)

            ## Links to TAIR & PARE gene information, plot location
            ## TAIR
            tair.link <- paste0(tair.prefix, 
                                pydeg_cat$ID,">",
                                pydeg_cat$ID, "</a>")

            ## PARE link
            pare.link <- paste0(pare.prefix, 
                                pydeg_cat$ID,
                                pare.suffix,
                                fa("up-right-from-square", width=icon_size, height=icon_size,
                                   fill = "steelblue"),
                                "</a>")

            ## Plot links
            G.dplot.link <- paste0(plot.prefix,'"',
                                   gene_plot_location,'"',
                                   plot.suffix1, "'",
                                   gene_plot_location, "'",
                                   plot.suffix_widthG,'">',
                                   fa("chart-area", width=icon_size, height=icon_size,
                                    fill = "steelblue"),
                                 "</a>")

            P.dplot.link <- paste0(plot.prefix,'"',
                                   peak_plot_location,'"',
                                   plot.suffix1, "'",
                                   peak_plot_location, "'",
                                   plot.suffix_widthP,'">',
                                   fa("chart-column", width=icon_size, height=icon_size,
                                    fill = "steelblue"),
                                 "</a>")

            ## Add links to data table
            pydeg_cat$PARE <- pare.link
            pydeg_cat$TAIR <- tair.link
            pydeg_cat$Gene_plot <- G.dplot.link
            pydeg_cat$Peak_plot <- P.dplot.link

            ## Convert discrete columns to factors
            pydeg_cat$category_1 <- factor(pydeg_cat$category_1, levels= c(1:4,0))
            pydeg_cat$category_2 <- factor(pydeg_cat$category_2, levels= c(LETTERS[1:3]))
            pydeg_cat[,":="(MorePeaks=ifelse(MorePeaks==0,"No","Yes"))]
            pydeg_cat[,":="(rep_gene=ifelse(rep_gene==0,"No","Yes"))]
            pydeg_cat$comparison <- i.comp

            itag <- paste(i.comparison,"-Cat_",cat_label)
            pydeg_dt_list[[itag]] <- pydeg_cat
        }
    }
  }
}
pydeg_dt <- do.call(rbind,pydeg_dt_list)

## Fix comparison as factor
pydeg_dt$comparison <- factor(as.character(pydeg_dt$comparison))

## Combine categories 1 and 2 into a single col
pydeg_dt$cat1_2 <- paste(pydeg_dt$category_1,pydeg_dt$category_2,sep = "-")

## Reorder cols
pydeg_dt <- dplyr::select(pydeg_dt,
              all_of(c(i.cols.sort[1:4],"TAIR",i.cols.sort[5:10],"PARE",
                       i.cols.sort[11:length(i.cols.sort)],"Index")))

## Change column name for html report
setnames(pydeg_dt,i.cols.sort, i.cols.export)					   

## ---- summary ----
ifile <- "Peak_counts_BeforeAfter_Filtering"
sumFile <- file.path(summary_dir,ifile)
sum.df <- read.table(sumFile, header=TRUE)

sum.df$Comparison <- as.factor(sum.df$Comparison)
i.comparisons <- NULL
mfs <- NULL
confs <- NULL
for (i in seq_len(nrow(sum.df))) {
    ## Format comparison
    i.comparison <- sum.df$Comparison[i]
    i.samples <- unlist(strsplit(gsub("_and_","-",i.comparison),
                                 split="-"))
    i.samples2 <- NULL
    for (j in i.samples) {
        i.sel <-  sample_list %in% j
        names(j) <- names(sample_list)[i.sel]
        i.samples2 <- c(i.samples2,j)
        i.comp <- mergeVector(i.samples2)
    }
    i.comp <- mergeVector(ivec=i.samples2, sep2 = " AND <br/>", single=TRUE)
    i.comparisons <- c(i.comparisons,i.comp)

    ## Split settings
    i.settings <- unlist(strsplit(sum.df$Settings[i],"-"))
    multiplicative_factor <- i.settings[1]
    conf_level <- i.settings[2]
    mfs <- c(mfs, multiplicative_factor)
    confs <-  c(confs, conf_level)
}

## Add new columns
sum.df$MF <- mfs
sum.df$conf <- confs
sum.df$Comparison <- i.comparisons

## Order and remove column, 'Settings'
sum.df2 <- sum.df[,c(6,8,9,1:5)]

##Change column names for html report
colnames(sum.df2) <- c("Comparison",
                       "Mult.<br/>Factor",
                      "Conf.<br/>level",
                      "Total peaks<br/>pyDegradome",
                      "No peaks<br/>Filtered",
                      "No peaks<br/>Class. 1<br/> Cat. 1",
                      "No peaks<br/>Class. 2<br/> Cat. A",
                      "No peaks<br/>Cat. 1 &<br/> Cat. A")

sum.dt <- DT::datatable(sum.df2,
                        escape = FALSE,
                      rownames = FALSE,
                      width="100%",
                      ## height = "10px",
                      filter = "top",
                      extensions = c("FixedColumns","Buttons","RowGroup"),
                      options = list(
                          scrollX = TRUE,
                          fixedColumns = list(leftColumns = 1, rightColumns = 1),
                          dom = 'Bfrtip',
                          buttons = c('copy', 'csv', 'excel'),
                          pageLength = 150
                      )) %>%
    formatStyle(
        colnames(sum.df2[ncol(sum.df2)]),
        backgroundColor = '#f080801a',
        fontWeight = 'bold'
    ) %>%
    formatStyle(
        colnames(sum.df2[ncol(sum.df2)-2]),
        backgroundColor = '#f080801a'
    ) %>%
    formatStyle(
        colnames(sum.df2[ncol(sum.df2)-1]),
        backgroundColor = '#add8e61a'
    ) %>%
    formatStyle(
        colnames(sum.df2[ncol(sum.df2)-3]),
        backgroundColor = '#add8e61a',
        fontWeight = 'bold'
    ) %>%
    formatStyle(
        colnames(sum.df2[1]),
        fontWeight = 'bold'
    )
sum.dt

## ---- general-classification ----
#Create DT object
i.dt <- DT::datatable(dplyr::select(pydeg_dt,all_of(cols.gen)),
                      escape = FALSE,
                      rownames = FALSE,
                      width="100%",
                      filter = "top",
                      extensions = c("FixedColumns","Buttons","RowGroup"),
                      options = list(
                          scrollX = TRUE,
                          scrollY = "500px",
                          scroller = TRUE,
                          fixedColumns = list(leftColumns = 4, rightColumns = 1),
                          columnDefs = list(list(visible=FALSE, targets="Comparison")),
                          rowGroup = list(dataSrc = length(cols.gen)-1),
                          dom = 'Bfrtip',
                          buttons = c('copy', 'csv', 'excel'),
                          pageLength = 150
                      )
                      ) %>%
    formatStyle(
    c("Peaks 2>"),
    color = styleEqual(
      c("No","Yes"), c("black", "Crimson")
    )) %>%
  formatStyle(
    "Feature",
    background = styleEqual(
      c("3UTR","CDS","5UTR"), c("#90ee901a", "#f080801a", "#add8e61a")
    )
  ) %>%
    formatRound(cols2round[5], 2)
i.dt

## ---- classification-metrics ----
i.dt2 <- DT::datatable(dplyr::select(pydeg_dt,all_of(cols.metrics)),
                       escape = FALSE,
                       rownames = FALSE,
                       width="100%",
                       filter = "top",
                       extensions = c("FixedColumns","Buttons","RowGroup"),
                       options = list(
                           scrollX = TRUE,
                           scrollY = "500px",
                           scroller = TRUE,
                           fixedColumns = list(leftColumns = 3),
                           columnDefs = list(list(visible=FALSE, targets="Comparison")),
                           rowGroup = list(dataSrc = length(cols.metrics)-1),
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 150)
                       ) %>%
    formatRound(cols2round[1:4], 2)
i.dt2

## ---- coordinates ----
i.dt3 <- DT::datatable(dplyr::select(pydeg_dt,all_of(cols.coordinates)),
                       escape = FALSE,
                       rownames = FALSE,
                       width="100%",
                       filter = "top",
                       extensions = c("FixedColumns","Buttons","RowGroup"),
                       options = list(
                           scrollX = TRUE,
                           scrollY = "500px",
                           scroller = TRUE,
                           fixedColumns = list(leftColumns = 3),
                           columnDefs = list(list(visible=FALSE, targets="Comparison")),
                           rowGroup = list(dataSrc = length(cols.coordinates)-1),
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 150)
                       )
i.dt3

## ---- summary-miRtargets ----
ifileA <- paste0("miRNA_seq/output/miRNA_targets_MF-", i.MF, "_iConf-", i.conf_f, ".txt")
mirFileA <- file.path(supp_data_dir, ifileA)

if (file.exists(mirFileA)) {
    miRtargetsA <- tryCatch(read.table(mirFileA, comment.char ="#", header=FALSE),
                           error = function(e) {
                               message(#
                                  cat("Peak sequences did not aligned with known miRNA species"))
                               message(e)
                           })
    
    if (class(miRtargetsA)=="data.frame") {
        miRtargetsA <- within(miRtargetsA, {
            ID <- ave(V1, V1, FUN=seq_along)
        })
        miRtargetsWideA <- reshape(miRtargetsA, direction  = "wide", idvar="ID", timevar="V1")
        colnames(miRtargetsWideA) <- mgsub::mgsub(colnames(miRtargetsWideA),c("V2.",":"),c("",""))
        miRtargetsWideA <- miRtargetsWideA[with(miRtargetsWideA,
                                              order(Comparison,
                                                    Transcript)),]
        ## Change format of score column
        miRtargetsWideA$Score <- as.numeric(miRtargetsWideA$Score)
        
        ## Format columns and add alignment image
        miRtargetsWideA <- dplyr::select(miRtargetsWideA,
                                         all_of(c("Transcript","miRNA",
                                                  "Score","Index","Comparison")))
        attach(miRtargetsWideA)
        img_file <- file.path("Alignment/global",Comparison,
                              paste0(paste("Aln_global",gsub("\\.","_",Transcript),
                                    miRNA, i.conf_f, i.MF, sep = "_"),".png"))
        detach(miRtargetsWideA)
        alignment.link <- paste0(plot.prefix,'"',
                                 img_file,'"',
                                 plot.suffix1, "'",
                                 img_file, "'",
                                 plot.suffix_widths,'">',
                                 fa("align-justify", width=icon_size, height=icon_size,
                                    fill = "steelblue"), "</a>")
        miRtargetsWideA$Alignment_global <- alignment.link
    }
}

ifileB <- paste0("miRNA_seq/output/mirmap_miRNA_targets_MF-", i.MF, "_iConf-", i.conf_f, ".txt")
mirfileB <- file.path(supp_data_dir, ifileB)

if (file.exists(mirfileB)) {
    miRtargetsB <- tryCatch(read.table(mirfileB, comment.char ="#", header=FALSE),
                           error = function(e) {
                               message(#
                                  cat("Peak sequences did not aligned with known miRNA species"))
                               message(e)
                           })
    
    if (class(miRtargetsB)=="data.frame") {
        miRtargetsB <- within(miRtargetsB, {
            ID <- ave(V1, V1, FUN=seq_along)
        })
        miRtargetsWideB <- reshape(miRtargetsB, direction  = "wide", idvar="ID", timevar="V1")
        colnames(miRtargetsWideB) <- mgsub::mgsub(colnames(miRtargetsWideB),c("V2.",":"),c("",""))
        miRtargetsWideB <- miRtargetsWideB[with(miRtargetsWideB,
                                              order(Comparison,
                                                    Transcript)),]

        ## Change format of score column
        miRtargetsWideB$Score <- as.numeric(miRtargetsWideB$Score)
        
        ## Format columns and add alignment image
        miRtargetsWideB <- dplyr::select(miRtargetsWideB,
                                         all_of(c("Transcript","miRNA","Score",
                                                  "Index","Comparison")))
        attach(miRtargetsWideB)
        img_file <- file.path("Alignment/mirmap",Comparison,
                              paste0(paste("Aln_mirmap",gsub("\\.","_",Transcript),
                                    miRNA, i.conf_f, i.MF, sep = "_"),".png"))
        detach(miRtargetsWideB)
        mirmap.link <- paste0(plot.prefix,'"',
                              img_file,'"',
                              plot.suffix1, "'",
                              img_file, "'",
                              plot.suffix_widthS,'">',
                              fa("align-left", width=icon_size, height=icon_size,
                                 fill = "steelblue"), "</a>")
        miRtargetsWideB$Alignment_mirmap <- mirmap.link
    }
}

## Merge dataframes
if (exists("miRtargetsWideA") && exists("miRtargetsWideB")) {
    attach(miRtargetsWideA)
    miRtargetsWideA$indx <- paste0(Transcript, miRNA, Comparison)
    detach(miRtargetsWideA)

    attach(miRtargetsWideB)
    miRtargetsWideB$indx <- paste0(Transcript, miRNA, Comparison)
    detach(miRtargetsWideB)

    miRtargetsWide <- merge(miRtargetsWideA,
                            miRtargetsWideB,
                        by = "indx", all = TRUE)[,-1]

    ## Merge common columns
    cols <- c("Transcript", "miRNA", "Index", "Comparison")
    for (col in cols) {
        isel <- grepl(col, colnames(miRtargetsWide))
        df_tmp <-  miRtargetsWide[,isel]
        icol <- NULL
        for (i in seq_len(nrow(df_tmp))) {
            tmp <- as.character(df_tmp[i,])
            tmp <- tmp[!is.na(tmp)]
            if (length(tmp)==1) {
                icol <- c(icol,tmp)
            } else if (tmp[1]==tmp[2]) {
                icol <- c(icol,tmp[1])
            } else {
                icol <- NA
            }
        }
        miRtargetsWide <- miRtargetsWide[,!isel]
        miRtargetsWide$tmp <- icol
        setnames(miRtargetsWide, "tmp", col)
    }
    
} else if (exists("miRtargetsWideA")) {
    miRtargetsWide <- miRtargetsWideA
} else {
    miRtargetsWide <- miRtargetsWideB
}

## Format comparison column
i.comparisons <- NULL
for (i in seq_len(nrow(miRtargetsWide))) {
    ## Format comparison
    i.comparison <- miRtargetsWide$Comparison[i]
    i.samples <- unlist(strsplit(gsub("_and_","-",i.comparison),
                                 split="-"))
    i.samples2 <- NULL
    for (j in i.samples) {
        i.sel <-  sample_list %in% j
        names(j) <- names(sample_list)[i.sel]
        i.samples2 <- c(i.samples2,j)
        i.comp <- mergeVector(i.samples2)
    }
    i.comp <- mergeVector(ivec=i.samples2, sep2 = " AND <br/>", single=TRUE)
    i.comparisons <- c(i.comparisons,i.comp)
}
miRtargetsWide$Comparison <- i.comparisons

## Merge link to peak plot
cols_sub <- c("Transcript",
              "Gene\nplot",
              "Peak\nplot",
              "Peak(T):Max(C)",
              "Feature",
              "Gene name",
              "Description", 
              "Comparison",
              "Index")
pydeg_dt_sub <- dplyr::select(pydeg_dt, all_of(cols_sub))
miRtargetsWide <- merge(dplyr::select(miRtargetsWide,
                                      c("miRNA", "Alignment_global",
                                        "Score.x",
                                        "Alignment_mirmap",
                                        "Score.y",
                                        "Index")),
                        pydeg_dt_sub,
                        by = "Index", all.x = TRUE)[,-1]

## Select columns for output
miRtargetsWide <- dplyr::select(miRtargetsWide,c("Transcript",
                                                 "miRNA",
                                                 "Alignment_global",
                                                 "Score.x",
                                                 "Alignment_mirmap",
                                                 "Score.y",
                                                 cols_sub[2:8]))
## Add missing-alignment-icon
no_link <- paste0("<a  target=_blank >",
             fa("ban", width=icon_size, height=icon_size,
                fill = "grey"), "</a>")

isel <- is.na(miRtargetsWide$Alignment_global)
miRtargetsWide$Alignment_global[isel] <- no_link

isel <- is.na(miRtargetsWide$Alignment_mirmap)
miRtargetsWide$Alignment_mirmap[isel] <- no_link

## Change alignment column name
aln_cols <- c("Alignment_global","Alignment_mirmap")
setnames(miRtargetsWide,
         aln_cols,
         gsub("_","\n",aln_cols))	


miRtargetsWide <- miRtargetsWide[with(miRtargetsWide, order(Comparison, -Score.x, Score.y)),]

## Change name of score columns
setnames(miRtargetsWide, "Score.x", "Alignment\nscore")
setnames(miRtargetsWide, "Score.y", "DeltaG\nopen")
                                 
miR.dt <- DT::datatable(miRtargetsWide,
                        escape = FALSE,
                        rownames = FALSE,
                        width="100%",
                        filter = "top",
                        extensions = c("FixedColumns","Buttons","RowGroup"),
                        options = list(
                            scrollX = TRUE,
                            scrollY = "500px",
                            scroller = TRUE,
                            fixedColumns = list(leftColumns = 1, rightColumns = 1),
                            columnDefs = list(list(visible=FALSE, targets="Comparison")),
                            rowGroup = list(dataSrc = ncol(miRtargetsWide)-1),
                            dom = 'Brtip',
                            buttons = c('copy', 'csv', 'excel'),
                            pageLength = 150
                        )) %>%
    formatStyle(
        "Feature",
        background = styleEqual(
            c("3UTR","CDS","5UTR"), c("#90ee901a", "#f080801a", "#add8e61a")
        )
    ) %>%
    formatRound(cols2round[5:6], 2)
miR.dt 
