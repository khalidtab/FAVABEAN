suppressMessages(library("optparse"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Path to Robject containing the error model and path to the fastq files", metavar="Path to Robject containing the error model and path to the fastq files"),
  make_option(c("-r", "--rdirection"), type="character", default=NULL, help="R1 or R2", metavar="R1 or R2"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="The number of cores to use in this script", metavar="The number of cores to use in this script")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

load(opt$input)

cores = as.numeric(opt$cores)
rDirection = opt$rdirection

# Denoising the sequences
dadaDenoised = dada(theFiles, err=DADA2errors, multithread=cores, pool=TRUE, verbose=2)
print("Done with denoising. Will write file now…")
rm(opt,cores)

save(dadaDenoised, theFiles, inputFolder,file = paste0(inputFolder,"_DADA2Denoise-",rDirection,".RData"), envir = .GlobalEnv)
