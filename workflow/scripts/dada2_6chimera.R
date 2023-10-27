suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input RObject file", metavar="input RObject file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="output file"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="cores", metavar="cores")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



inputfile = opt$input
outputFile = opt$output
cores = as.numeric(opt$cores)

load(inputfile)

seqtab_noChim =  removeBimeraDenovo(seqtab_combined, method = "pooled", multithread =  cores, verbose=TRUE)

print("Chimeras removed")


save(list = list(seqtab_noChim), file = outputFile, envir = .GlobalEnv)
