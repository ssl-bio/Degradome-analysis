##Description:  Wrapper for running Rmd files

## Libraries
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(PostPyDeg))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fontawesome))

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
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    i.conf_f <- gsub("\\.","_",i.conf)
    my.input <- here("Scripts/R/05-Report.Rmd")
    out.name <- paste0("Top candidate Targets from Degradome Analysis MF-",i.MF,
                       "_Conf-",i.conf_f, "_",ibase,
                       ".html")
    if(!file.exists(file.path(report_dir,out.name))) {
        rmarkdown::render(my.input,
                          params = list(i.MF=i.MF,
                                        i.conf=i.conf,
                                        output.dir=report_dir),
                          output_file = file.path(report_dir,out.name))
    }
}
