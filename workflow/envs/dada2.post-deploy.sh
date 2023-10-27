conda activate $CONDA_PREFIX
Rscript -e 'install.packages("BiocManager",repos = "http://cran.us.r-project.org"); BiocManager::install("dada2"); install.packages(c("readr","magrittr","optparse","tidyr","jsonlite","furrr","future"),repos="http://cran.us.r-project.org")'
