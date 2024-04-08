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
forwardfiles = inputfiles %>% grepl(myparameter,x=.) %>% inputfiles[.] %>% grepl("seqtab_forwardOnly",x=.)  %>% inputfiles[.]
pairedfiles = inputfiles[-which(inputfiles %in% forwardfiles)]

process_file <- function(path) {
  newTable <- read_tsv(path) %>% 
    as.data.frame() 
  
  rownames(newTable) <- newTable[, 1]
  newTable[, 1] <- NULL
  
  newTable = as.matrix(newTable)
  return(newTable)
}

forwardfiles = base::sapply(forwardfiles, process_file,simplify=TRUE)
pairedfiles  = base::sapply(pairedfiles , process_file,simplify=TRUE)

print("Merging sequence tables, then performing chimera removal")

forwardfiles = forwardfiles %>% mergeSequenceTables(tables=.,tryRC = TRUE) %>% removeBimeraDenovo(., multithread =  TRUE, verbose=TRUE, method = "pooled")
pairedfiles  = pairedfiles  %>% try(mergeSequenceTables(tables=.,tryRC = TRUE),silent=FALSE) %>% try(removeBimeraDenovo(., multithread =  TRUE, verbose=TRUE, method = "pooled"))

print("Chimeras removed")

save(forwardfiles, pairedfiles, file = paste0("data/favabean/",myparameter,"_chimeraRemoved.RObjects"), envir = .GlobalEnv)
