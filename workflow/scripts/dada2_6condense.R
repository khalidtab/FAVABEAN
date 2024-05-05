suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file after chimera removal", metavar="input file after chimera removal"),
  make_option(c("-r", "--region"), type="character", default=NULL, help="input region", metavar="input region"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="output file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



inputFile = opt$input
cores = opt$cores

load(inputFile)


print("starting collapse of ASVs…")
condensed_table = collapseNoMismatch(inputfile4, minOverlap = 20, identicalOnly = FALSE, band = -1, verbose = TRUE)

SampleIDs = rownames(condensed_table)
condensed_table = as.data.frame(condensed_table)
condensed_table = cbind(SampleIDs,condensed_table)
condensed_table = t(condensed_table) %>% as.data.frame(.)
condensed_table = cbind(rownames(condensed_table),condensed_table)

write_tsv(condensed_table,opt$output,col_names = FALSE)
