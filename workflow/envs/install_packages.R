#!usr/bin/env Rscript


#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2")

BiocManager::install("optparse",ask=FALSE,force=TRUE,update=FALSE)
BiocManager::install(c("readr", "magrittr", "tidyr", "jsonlite","data.table","dplyr"),ask=FALSE,force=TRUE,update=FALSE)
BiocManager::install(c("mgcv", "DelayedArray", "SummarizedExperiment", "GenomicAlignments", "ShortRead", "ggplot2", "dada2"),ask=FALSE,force=TRUE,update=FALSE)
