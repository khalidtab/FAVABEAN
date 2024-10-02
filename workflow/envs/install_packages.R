#!usr/bin/env Rscript


#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2")

BiocManager::install("optparse",ask=FALSE,force=TRUE,update=FALSE)
BiocManager::install(c("readr", "magrittr", "tidyr", "jsonlite","data.table"),ask=FALSE,force=TRUE,update=FALSE)
BiocManager::install("dada2",,ask=FALSE,force=TRUE,update=FALSE)
