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
for (i in seq_along(MF_list)) {
    i.MF <- MF_list[i]
    i.conf <- conf_list[i]
    i.conf_f <- gsub("\\.","_",i.conf)
    my.input <- here(paste0("Scripts/R/05-Report_",ibase,".Rmd"))
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
