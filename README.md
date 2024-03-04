# Degradome analysis


## Description

Pipeline to map mRNA degradome sequences, compare two conditions and obtain a list of transcripts containing regions with significant differences. All the steps are grouped into two scripts, 01-Degradome.sh and 02-Degradome.sh. The former's main job is to map the sequencing data and call for a third party script (`PyDegradome` <a href="#citeproc_bib_item_1">Gaglia, Rycroft, and Glaunsinger 2015</a>) which will compare two samples in order to identify regions with significant differences (*peaks*). The latter script will then, associate these regions with known transcripts, and classify them using ratios of read numbers in the peak region and elsewhere (within and between samples).


## Environment setup

A number of different programs are used and these are installed in three miniconda environments:

-   **pydeg\_map**: Used for most of the steps of the first script, `01-Degradome.sh`
-   **pydeg\_run**: Used mainly to run the PyDegradome script.
-   **pydeg\_R**: Used to run all the R scripts.
-   **pydeg\_aux**: Used to install packages that cause conflicts: `salmon`

The following commands will create all the anaconda environments mentioned above and install the required packages within them. Linux dependencies needed to install and setup R are indicated at the beginning.

```sh
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

conda install -c bioconda bioawk fastqc trim-galore bedtools seqkit bowtie2 picard samtools biopython qualimap htseq deeptools entrez-direct sra-tools -y
conda install -c bioconda agat #Should be done separately to avoid hanging
conda install -c anaconda pandas pilow tk pip -y
conda install -c conda-forge dash dash-bootstrap-components pdf2image pigz -y
conda install -c plotly plotly plotly_express -y
pip install dash-dangerously-set-inner-html dash_bootstrap_templates nbib dash_breakpoints

set +u
conda deactivate

# Create R environment
conda create --name pydeg_R python=3.11.4 -y
source activate pydeg_R
set -u

conda install -c conda-forge pkg-config=0.29.2 r-curl=5.0.1 r-base=4.3.1 -y

set +u
conda deactivate

# Create environment for running pydegradome
conda create --name pydeg_run python=2 -y
source activate pydeg_run
set -u

conda install -c bioconda biopython tophat viennarna unittest2 dendropy -y
conda install -c anaconda scipy pandas pillow tk -y

set +u
conda deactivate

# Create auxiliary environment
conda create --name pydeg_aux python=3 -y
source activate pydeg_aux
set -u

conda install -c bioconda salmon=1.10.3 -y
```


### R setup

The following commands should be run on the environment pydeg\_R. They will install bioconductor and all the packages required for the analysis of the data. Note that the version of bioconductor should be compatible with that of R (See [bioconductor installation instructions](https://www.bioconductor.org/install/))

```R
## Bioconductor install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

## Install all packages
ipackages <- c("Biostrings", "biomaRt", "bookdown", "Cairo",
               "ChIPpeakAnno", "DESeq2", "devtools", "DT",
               "GenomicFeatures", "Gviz", "RColorBrewer", "RVenn",
               "biomaRt", "data.table", "doParallel", "dplyr",
               "ensembldb", "filesstrings", "fontawesome", "ggplot2",
               "gsubfn", "here", "knitr", "magrittr", "mgsub", "optparse", "pbapply", "purrr", "reshape2", "rmarkdown", "rtracklayer", "seqinr", "stringr", "tidyverse")
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


## Running the plotly dash app

Some modifications were included in order to summarize the results of the analysis in a plotly dash app. The second script (`02-Degradome.sh`) now produces a new folder `04-Dash_app`, with all the files needed to display a report in a web browser. To run the app is necessary to activate the environment, `pydeg_map`, and from within the folder run the following command:

```bash
python app.py
```

Then, the report can be analyzed on a web browser on the address, `localhost:8050`

A live demo can be found at <https://sslbio.pythonanywhere.com/>


## Modifications


### Regarding the files [August - October 2023]

To produce a plotly dash app the following files were included

-   `Env_variables/PostPydeg_factor_description.tsv`: Tab-separated file with details of the analysis.
-   `Env_variables/custom.css`: css file specifying some style options for the report
-   `Scripts/sh_py/06-Process_data.py`: Script for processing the data required for the report
-   `Scripts/sh_py/app.py`: Plotly dash app.

The files above are in an early stage and some work is needed to make the process more efficient, particularly.

-   [ ] Improve the generation of the `Env_variables/PostPydeg_factor_description.tsv` file. Currently the description was written in emacs org-mode, exported to html format and copied to the `description` column in the `tsv` file.
-   [ ] Improve the generation of files and variables (06-Process\_data.py). Currently peak plots generated by the script, `03-Drawing_Dplots.R` are converted from `pdf` to `jpg` duplicating the number of plots.
    -   One option is to delete the `pdf` files and map the links in the reports generated in R to the location of the `jpg` equivalents
-   [X] Improve the quality of the report (app.py)


### Fixing the plotly dash app [March 2024]

A previous version of the plotly dash app would fail to run in the following cases

-   A project name different from `Zhang-2021` was used
-   One or both miRNA alignments were not produced The former was caused because one of the app files (*i.e.* `test_cases.py` ) had the project name hard coded. Now this is replaced by the project name used when running `06-Process_data.py`. The latter was solved by including a number of tests for missing alignment data. The missing data was observed when running the analysis on a sub-sample of the original data which caused the top categories to be absent.


## Further details

For further details see this [post](https://ssl-blog.netlify.app/posts/degradome-analysis/degradome-code/)


## References
  <div class="csl-entry"><a id="citeproc_bib_item_1"></a>Gaglia, Marta Maria, Chris H. Rycroft, and Britt A. Glaunsinger. 2015. “Transcriptome-Wide Cleavage Site Mapping on Cellular mRNAs Reveals Features Underlying Sequence-Specific Cleavage by the Viral Ribonuclease SOX.” Edited by Pinghui Feng. <i>Plos Pathogens</i> 11 (12): e1005305. doi:<a href="https://doi.org/10.1371/journal.ppat.1005305">10.1371/journal.ppat.1005305</a>.</div>
</div>
