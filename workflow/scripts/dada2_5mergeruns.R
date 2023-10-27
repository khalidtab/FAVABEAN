suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("future"))
suppressMessages(library("furrr"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input folder with all the seqtab files", metavar="input folder with all the seqtab files"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file", metavar="output with all the seqtab files"),
  make_option(c("-r", "--region"), type="character", default=NULL, 
              help="region of interest", metavar="region of interest")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

inputFolder = opt$input
outputFile = opt$output
region = opt$region

# Set up parallel processing
plan(multisession)

# List all files with "seqtab.tsv" pattern and filter for those that contain the 'region' variable
seqtabs <- list.files(inputFolder, pattern="seqtab.tsv", full.names = TRUE, recursive = TRUE) %>% 
  grep(pattern = region, value = TRUE)

# Define a function to read and reshape each table
read_and_reshape <- function(x) {
  cat("Currently merging table", x, "\n")
  read_tsv(x) %>%
    mutate(across(-V1, as.character), source = x) %>% # Convert all non-V1 columns to character and add source column for the filename
    pivot_longer(cols = -c(source, V1), names_to = "sequence", values_to = "abundance")
}

# Parallelize the reading and reshaping operation
combinedTables <- future_map_dfr(seqtabs, read_and_reshape)

# Pivot operations: longer then wider
seqtab_combined <- combinedTables %>%
  pivot_wider(names_from = "sequence", values_from = "abundance", values_fill = "0") %>%
  as.data.frame() %>%
  mutate(source = make.unique(as.character(source))) %>%
  { rownames(.) <- .$source; . } %>%
  select(-source)

rm(combinedTables, read_and_reshape, seqtabs)

print("All tables combined")

save(list = ls(all.names = TRUE), file = paste0(inputFolder, region, "_mergedRuns.RObjects"), envir = .GlobalEnv)
