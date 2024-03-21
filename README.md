# Overview
The code provided in this repository allows the computation of some main statistics and the production of the figures included in the article "Idiosyncratic and generic single nuclei and spatial transcriptional patterns in papillary and anaplastic thyroid cancers".

It starts with the pre-processed R objects available on Zenodo under DOI 10.5281/zenodo.842152.

# System requirements
## OS requirements
The containers were created and the analyses were run on Ubuntu 20.04 LTS. Systems that can run apptainer are recommended. Alternatively, systems that can run R 4.1.0 will suffice. For Apptainer compatible systems, see https://apptainer.org/docs/admin/main/installation.html#installation-on-windows-or-mac
## Software dependencies
This code was run in R4.1.0, with the list of packages provided below installed.

The R session was running in apptainer version 1.0.2.
## R packages required
* Seurat
* ggplot2
* catecolors
* harmony
* rentrez
* org.Hs.eg.db
* ComplexHeatmap
* circlize
* stringr
* CellChat
* fgsea
* msigdbr


# Installation guide
## Apptainer installation
The Apptainer container can be created on any system where Apptainer can be installed. To install apptainer, see : https://apptainer.org/docs/admin/main/installation.html
## Container creation
The container can be created by downloading the apptainer_creation.def file and using the command the command below in a linux terminal with apptainer installed:

sudo apptainer build apptainer_igp.sif apptainer_creation.def

## Alternative installation
To run this R code without Apptainer, installing R 4.1.0 and the packages listed above will suffice. See https://cran.r-project.org/bin/
## Installation time
Installation time can vary depending on the system. It typically lasts over an hour.

# Demo
To run the R scripts provided in this repository, start the container, change the working directory to the directory containing the scripts and load the configuration scripts:

apptainer exec apptainer_igp.sif R

setwd("/mnt/iribhm/ngs/ST/article/R/")

source("setup.R")

source("constants.R")

source("utils.R")

source("graphics.R")

Starting the container and sourcing the configuration scripts will prepare the R session to run any scripts contained in this repository. It should last less than one minute. 

# Usage instruction
Any of the R scripts included in this repository can be run from any installation of R4.1.0 containing the packages mentioned above. Once your R session is started, the scripts are designed to work from the directory containing the script files (i.e. this repository), sitting next to the folders containing the inputs and outputs. Those input and output folders and their subfolders structure are listed in the constants.R script under the "DIRECTORIES" heading. Running all scripts from this repository should reproduce the figures and tables from the article.

The time required to run all scripts included here vary widely depending on the system resources. It will last several hours minimum.
