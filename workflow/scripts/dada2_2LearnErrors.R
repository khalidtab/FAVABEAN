# Have the same script to learn errors of forward and reverse, but it depends on the input

suppressMessages(library("optparse"))
suppressMessages(library("dada2"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Folder with the trimmed and filtered FASTQ files from DADA2", metavar="Folder with the trimmed and filtered FASTQ files from DADA2"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="Number of cores", metavar="Number of cores"),
  make_option(c("-r", "--rdirection"), type="character", default=NULL, help="Is this to be run on R1 or R2 reads", metavar="Is this to be run on R1 or R2 reads")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

inputFolder = opt$input
rDirection = opt$rdirection
cores = as.numeric(opt$cores)

theFiles = sort(list.files(paste0(inputFolder,"/filtered"), pattern=paste0(rDirection,".*fastq.*"), full.names = TRUE, recursive = TRUE))

print("DADA2 Filtered and trimmed fastq files found")
print(theFiles)

print(paste("Will now run DADA2 Error learning. Number of cores to be utilized:",cores))

DADA2errors = learnErrors(theFiles, randomize=TRUE, multithread = cores, verbose = 2)

print("DADA2 Errors learned")

rm(cores,opt)

save(DADA2errors, inputFolder, theFiles, file = paste0(inputFolder,"_DADA2Errors-",rDirection,".RData"), envir = .GlobalEnv)

