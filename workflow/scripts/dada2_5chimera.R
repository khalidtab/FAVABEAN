suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input folder with the seqtab files to be combined", metavar="input folder with the seqtab files to be combined"),
  make_option(c("-p", "--parameter"), type="character", default=NULL, help="from all the seqfiles in that folder, any selection parameter you want to apply, eg, indication of the sequenced region", metavar="indication of the sequenced region"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="output file"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="cores", metavar="cores")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}


outputFile = opt$output
cores = as.numeric(opt$cores)

inputfile = opt$input
myparameter = opt$parameter

inputfiles = list.files(inputfile,pattern = "seqtab", full.names = TRUE, recursive = TRUE)
inputfiles2 = inputfiles %>% grepl(myparameter,x=.) %>% inputfiles[.]

process_file <- function(path) {
  newTable = data.table::fread(path,sep="\t") %>% as.data.frame() # This is a faster way to read the files than read_tsv()
  
  rownames(newTable) <- newTable[, 1]
  newTable[, 1] <- NULL
  
  newTable = as.matrix(newTable)
  return(newTable)
}

inputfiles2 = base::sapply(inputfiles2, process_file,simplify=FALSE)
print("Merging sequence tables…")

if (length(inputfiles2) == 1){
  inputfiles2 = inputfiles2[[1]] %>% as.matrix(.)
} else {
inputfiles2 = inputfiles2 %>% mergeSequenceTables(tables=.,tryRC = TRUE)
}

print("Done with merging sequence tables. Table format modification to accomodate chimera removal starts now…")

seqtab_noChim =  removeBimeraDenovo(inputfiles2, multithread =  TRUE, verbose=TRUE, method = "consensus")

print("Chimeras removed")

save(seqtab_noChim, file = paste0("data/favabean/",myparameter,"_chimeraRemoved.RObjects"), envir = .GlobalEnv)