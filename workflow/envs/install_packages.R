#!usr/bin/env Rscript

BiocManager::install("optparse")
install.packages(c("readr", "magrittr", "optparse", "tidyr", "jsonlite","data.table"), repos="http://cran.us.r-project.org",verbose=TRUE,quiet=FALSE)
#BiocManager::install("dada2")

#install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2")