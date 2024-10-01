#!usr/bin/env Rscript


#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2")

BiocManager::install("optparse")
install.packages(c("readr", "magrittr", "optparse", "tidyr", "jsonlite","data.table"), repos="http://cran.us.r-project.org")
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("dada2")
