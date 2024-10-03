#!usr/bin/env Rscript

# Set a CRAN mirror to download packages from
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install core packages without updating other packages
BiocManager::install("optparse", ask=FALSE, force=TRUE, update=FALSE)
BiocManager::install(c("readr", "magrittr", "tidyr", "jsonlite", "data.table"), ask=FALSE, force=TRUE, update=FALSE)

# Install remotes to handle older versions
# Install the inflection package from CRAN
install.packages("inflection")

# Load the package into your R session
install.packages("remotes")

# Install the correct versions of MASS and Matrix
remotes::install_version("MASS", version = "7.3-55", repos = "https://cloud.r-project.org")
remotes::install_version("Matrix", version = "1.5-1", repos = "https://cloud.r-project.org")

# Install dada2 without updating other packages
BiocManager::install("dada2", ask=FALSE, force=TRUE, update=FALSE)
