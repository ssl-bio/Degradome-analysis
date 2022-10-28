## ---- setup ----
knitr::opts_chunk$set(echo = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      figure.width=14,
                      results='asis',
                      fig.pos = 'center') #out.width='\\textwidth'

# Set variables
i.conf_f <- gsub("\\.","_",params$i.conf)
i.MF <- params$i.MF

# Set variables
icon_size <- "25px"
pydeg_input_f <- file.path(pydeg_pooled_dir,
                           paste("Pooled",
                                 i.conf_f, "4", i.MF, sep="_"))
cols.gen <- c(i.cols.export[1:6],i.cols.export[24],"PARE",
              i.cols.export[7:11],"TAIR",i.cols.export[12:13])
cols.metrics <- c(i.cols.export[1:3],
                  i.cols.export[c(14:22,13)])
cols.coordinates <- c(i.cols.export[1:3],
                      i.cols.export[c(25:33,13)])								 
								 

                                        # Import data
pydeg_all <- fread(pydeg_input_f)
setDT(pydeg_all)
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
            
            ## Focus on the combination Cat1=1 and Cat2=A
            if(unique(pydeg_cat$category_1)=="1" &&
               unique(pydeg_cat$category_2)=="A") {
                Nplots <- length(pydeg_cat$category_1)
            } else {
                Nplots <- 5
            }

            ##Sort dataframe by txRatio
            ## setDF(pydeg__cat)
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

                                        # TAIR
            tair.link <- paste0(tair.prefix, 
                                pydeg_cat$ID,">",
                                pydeg_cat$ID, "</a>")

                                        # PARE link
            pare.link <- paste0(pare.prefix, 
                                pydeg_cat$ID,
                                pare.suffix,
                                fa("up-right-from-square", width=icon_size, height=icon_size,
                                   fill = "steelblue"),
                                "</a>")

                                        # Plot links
            G.dplot.link <- paste0(dplot.prefix,
                                 gene_plot_location,">",
                                 fa("chart-area", width=icon_size, height=icon_size,
                                    fill = "steelblue"),
                                 "</a>")

            P.dplot.link <- paste0(dplot.prefix,
                                 peak_plot_location,">",
                                 fa("chart-column", width=icon_size, height=icon_size,
                                    fill = "steelblue"), "</a>")

            ## Add links to data table
            pydeg_cat$PARE <- pare.link
            pydeg_cat$TAIR <- tair.link
            pydeg_cat$Gene_plot <- G.dplot.link
            pydeg_cat$Peak_plot <- P.dplot.link

                                        # Convert discrete columns to factors
            pydeg_cat$category_1 <- factor(pydeg_cat$category_1, levels= c(1:4,0))
            pydeg_cat$category_2 <- factor(pydeg_cat$category_2, levels= c(LETTERS[1:3]))
            pydeg_cat[,":="(miRNA=ifelse(miRNA==0,"No","Yes"))]
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

pydeg_dt <- dplyr::select(pydeg_dt,
              all_of(c(i.cols.sort[1:3],"TAIR",i.cols.sort[4:11],"PARE",
                       i.cols.sort[12:length(i.cols.sort)])))
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

## Order and remove columns
sum.df2 <- sum.df[,c(6,8,9,1:5)]

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
                          ## scrollY = "500px",
                          ## scroller = TRUE,
                          fixedColumns = list(leftColumns = 1, rightColumns = 1),
                          ## columnDefs = list(list(visible=FALSE, targets="Comparison")),
                          ## rowGroup = list(dataSrc = 7),
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
        ## fontSize = '75%',
        fontWeight = 'bold'
    )
sum.dt
## ---- general-classification ----
#Create DT object
i.dt <- DT::datatable(dplyr::select(pydeg_dt,all_of(cols.gen)),
                      escape = FALSE,
                      rownames = FALSE,
                      width="100%",
                      ## height = "10px",
                      filter = "top",
                      extensions = c("FixedColumns","Buttons","RowGroup"),
                      options = list(
                          scrollX = TRUE,
                          scrollY = "500px",
                          scroller = TRUE,
                          fixedColumns = list(leftColumns = 3, rightColumns = 1),
                          columnDefs = list(list(visible=FALSE, targets="Comparison")),
                          rowGroup = list(dataSrc = 15),
                          dom = 'Bfrtip',
                          buttons = c('copy', 'csv', 'excel'),
                          pageLength = 150
                      )
                      ) %>%
    formatStyle(
    c("miRNA","Peaks 2>"),
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

## ---- classification1-metrics ----

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
                           rowGroup = list(dataSrc = 12),
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 150)
                       ) %>%
    formatRound(cols2round[1:4], 2)
i.dt2

## ---- classification2-coordinates ----
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
                           rowGroup = list(dataSrc = 12),
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 150)
                       )
i.dt3

## ---- summary-miRtargets ----
ifile <- paste0("miRNA_targets_MF-", i.MF, "_iConf-", i.conf_f, ".txt")
mirFile <- file.path(summary_dir, ifile)
miRtargets <- read.table(mirFile, comment.char ="#", header=FALSE)
miRtargets <- within(miRtargets, {
  ID <- ave(V1, V1, FUN=seq_along)
})
miRtargetsWide <- reshape(miRtargets, direction  = "wide", idvar="ID", timevar="V1")
colnames(miRtargetsWide) <- c("ID", "Transcript", "Comparison", "miRNA")
miRtargetsWide <- miRtargetsWide[with(miRtargetsWide,
                                          order(Comparison,
                                                Transcript,
                                                miRNA)),]
miRtargetsWide <- as.data.frame(miRtargetsWide %>% group_by_at(vars(Transcript,Comparison)) %>%
  summarize_all(paste, collapse=","))

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

## Add new columns
miRtargetsWide$Comparison <- i.comparisons

## Merge link to peak plot
miRtargetsWide$indx <- paste0(miRtargetsWide$Transcript, miRtargetsWide$Comparison)
pydeg_dt_sub <- dplyr::select(pydeg_dt,c("Transcript","Peak\nplot","Comparison"))
pydeg_dt_sub$indx <- paste0(pydeg_dt_sub$Transcript, pydeg_dt_sub$Comparison)
miRtargetsWide <- merge(miRtargetsWide,
                        dplyr::select(pydeg_dt_sub,c("Peak\nplot","indx")),
                        by = "indx", all.x = TRUE)

## Select columns for output
miRtargetsWide <- dplyr::select(miRtargetsWide,c("Transcript",
                                                 "miRNA",
                                                 "Peak\nplot",
                                                 "Comparison"))

miRtargetsWide$miRNA <- gsub(",","<br/>",miRtargetsWide$miRNA)
miR.dt <- DT::datatable(miRtargetsWide,
                        escape = FALSE,
                      rownames = FALSE,
                      width="100%",
                      ## height = "10px",
                      ## filter = "top",
                      extensions = c("FixedColumns","Buttons","RowGroup"),
                      options = list(
                          scrollX = TRUE,
                          ## scrollY = "500px",
                          ## scroller = TRUE,
                          fixedColumns = list(leftColumns = 1, rightColumns = 1),
                          columnDefs = list(list(visible=FALSE, targets="Comparison")),
                          rowGroup = list(dataSrc = 3),
                          dom = 'Brtip',
                          buttons = c('copy', 'csv', 'excel'),
                          pageLength = 150
                      ))
miR.dt 
