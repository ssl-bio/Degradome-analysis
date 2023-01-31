# Degradome analysis


## Description

Pipeline to map mRNA degradome sequences, compare two conditions and obtain a list of transcripts containing regions with significant differences. All the steps are grouped into two scripts, 01-Degradome.sh and 02-Degradome.sh. The former's main job is to map the sequencing data and call for a third party script (PyDegradome <a href="#citeproc_bib_item_1">Gaglia, Rycroft, and Glaunsinger 2015</a>) which will compare two samples in order to identify regions with significant differences (*peaks*). The latter script will then, associate these regions with known transcripts, and classify them using ratios of read numbers in the peak region and elsewhere (within and between samples).


## Environment setup

A number of different programs are used and these are installed in three miniconda environments:

-   **pydeg\_map**: Used for most of the steps of the first script, 01-Degradome.sh
-   **pydeg\_run**: Used mainly to run the PyDegradome script.
-   **pydeg\_R**: Used to run all the R scripts.

The following commands will create all the anaconda environments mentioned above and install the required packages within them. Linux dependencies needed to install and setup R are indicated at the beginning.

```bash
# Install linux dependecies
sudo apt install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcurl4-openssl-dev libxml2-dev

eval "$(conda shell.bash hook)"

# Add additional channels
conda config --add channels conda-forge
conda config --add channels bioconda

# Create mapping environment
conda create --name pydeg_map python=3 -y
source activate pydeg_map

#Abort if get any error
set -eo pipefail

conda install -c bioconda bioawk fastqc trim-galore bedtools seqkit bowtie2 picard samtools biopython qualimap htseq deeptools salmon -y
conda install -c bioconda agat #Should be done separately to avoid hanging
conda install -c anaconda pandas pilow tk -y

set +u
conda deactivate

# Create R environment
conda create --name pydeg_R python=3 -y
source activate pydeg_R
set -u

conda install -c conda-forge pkg-config r-curl r-base -y

set +u
conda deactivate
# Create auxiliary environment
conda create --name pydeg_run python=2 -y
source activate pydeg_run

conda install -c bioconda tophat viennarna unittest2 dendropy -y
conda install -c anaconda biopython pandas pillow tk -y
```


### R setup

The following commands should be run on the environment pydeg\_R. They will install bioconductor and all the packages required for the analysis of the data. Note that the version of bioconductor should be compatible with that of R (See [bioconductor installation instructions](https://www.bioconductor.org/install/))

```r
## Bioconductor install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

## Install all packages
ipackages <- c("Biostrings", "bookdown", "Cairo", "ChIPpeakAnno", "DESeq2", "devtools", "DT",
	       "GenomicFeatures", "Gviz", "RColorBrewer", "RVenn", "biomaRt", "data.table",
	       "doParallel", "dplyr", "ensembldb", "filesstrings", "fontawesome", "ggplot2",
	       "gsubfn", "here", "knitr", "magrittr", "mgsub", "optparse", "pbapply", "purrr",
	       "reshape2", "rmarkdown", "rtracklayer", "seqinr", "stringr", "tidyverse")
BiocManager::install(ipackages)
```

1.  Custom package installation

    A number of commands used to process the data were wrapped in an R packaged named PostPydeg available at [github](https://github.com/ssl-bio/R_postpydeg.git). The commands to install it are described below.
    
    ```R
    library(devtools)
    devtools::document()
    devtools::build()
    devtools::install()
    ```


## File and Folder structure

The whole analysis pipeline uses and generates many files and folders. In an attempt to keep these in order, a specific file and folder structure was defined for both the input and the output. Most of these are defined in the variable setup file (contained in the `Env_variables` folder). Importantly, the output folder will have the following form, `Degradome-<project name>` and the sequencing files must be placed in the sub folder, `Degradome-<project name>/output_01/01-fastq` which should be created manually. As for the input folder, it could be assigned any name at the time of cloning. Execution of the analysis should be done at the root of the tree


## Running the analysis

In the commands below, note that two parameters are required, 'project name' and 'variable setup file'. The former should match the last part of the output folder (*e.g.* Degradome-Oliver-2022) and the latter can be named arbitrarily with the only requirement it must be located in the `Env_variables` folder.

```bash
# Arguments:
# /bin/bash 01-Degradome.sh <project name> <variable setup file>
# Example
/bin/bash 01-Degradome.sh Oliver-2022 Oliver-2022_vars.txt && \
    /bin/bash 02-Degradome.sh Oliver-2022 Oliver-2022_vars.txt
```


## Further details

For further details see this [post](https://ssl-blog.netlify.app/posts/degradome-analysis/degradome-code/)


## References
  <div class="csl-entry"><a id="citeproc_bib_item_1"></a>Gaglia, Marta Maria, Chris H. Rycroft, and Britt A. Glaunsinger. 2015. “Transcriptome-Wide Cleavage Site Mapping on Cellular mRNAs Reveals Features Underlying Sequence-Specific Cleavage by the Viral Ribonuclease SOX.” Edited by Pinghui Feng. <i>PLOS Pathogens</i> 11 (12): e1005305. doi:<a href="https://doi.org/10.1371/journal.ppat.1005305">10.1371/journal.ppat.1005305</a>.</div>
</div>
