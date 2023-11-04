suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input combined sequence file", metavar="input combined sequence file"),
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

seqtab_combined = read_tsv(inputfile,col_names= c("V1","sequence","abundance"), col_types = list("c","c","i"))

print(paste0("Starting chimera removal on file",inputfile))

seqtab_noChim =  removeBimeraDenovo(seqtab_combined, multithread =  TRUE, verbose=TRUE)

print("Chimeras removed")

write_tsv(seqtab_noChim,outputFile)
