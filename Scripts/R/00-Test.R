library(here)
library(stringr)
library(filesstrings)
## library(ggplot2)
library(biomaRt)
library(optparse)

## Parse arguments
option_list <- list(
    make_option(c("-d", "--wd"), type="character", default=NULL, 
                help="Working directory", metavar="character")
) 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Changing to directory", opt$wd)
