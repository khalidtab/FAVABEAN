suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file after chimera removal", metavar="input file after chimera removal"),
  make_option(c("-r", "--region"), type="character", default=NULL, help="input region", metavar="input region")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



inputFile = opt$input

seqtab_combined = read_tsv(inputFile,col_names= c("V1","sequence","abundance"), col_types = list("c","c","i"))

seqtab_combined = pivot_wider(seqtab_combined, names_from = "sequence", values_from = "abundance")
seqtab_combined = as.data.frame(seqtab_combined)
row.names(seqtab_combined) = seqtab_combined[,1]
seqtab_combined[,1] = NULL
seqtab_combined = as.matrix(seqtab_combined)

condensed_table = collapseNoMismatch(seqtab_combined, minOverlap = 20, orderBy = "abundance", identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)
SampleIDs = rownames(condensed_table)
condensed_table = as.data.frame(condensed_table)
condensed_table = cbind(SampleIDs,condensed_table)

region = opt$region

pathForCondensed_table = strsplit(inputFile,"favabean")[[1]][1]



write_tsv(condensed_table,
          paste0(pathForCondensed_table,"favabean/",region,"_condensed_ASVs.tsv"))
