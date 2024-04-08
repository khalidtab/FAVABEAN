suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Folder input with dada2 denoised files", metavar="Folder input with dada2 denoised files"),
  make_option(c("-f", "--fastq"), type="character", default=NULL, help="Folder input with dada2 filtered fastq files", metavar="Folder input with dada2 filtered fastq files")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

theFolder = opt$input
fastqFolder = opt$fastq

# Denoised files
dada2denoisedR1 = list.files(theFolder,pattern="_DADA2Denoise-R1.RData",full.names = TRUE)
dada2R1titles = basename(dada2denoisedR1) %>% strsplit(.,"_DADA2Denoise-R1.RData") %>% unlist(.) %>% strsplit(.,"dada2_") %>% unlist(.) %>% Filter(function(x) nchar(x) > 0, .)
dadaR1 = data.frame(title=dada2R1titles,denoiseR1=dada2denoisedR1)

dada2denoisedR2 = list.files(theFolder,pattern="_DADA2Denoise-R2.RData",full.names = TRUE)
dada2R2titles = basename(dada2denoisedR2) %>% strsplit(.,"_DADA2Denoise-R2.RData") %>% unlist(.) %>% strsplit(.,"dada2_") %>% unlist(.) %>% Filter(function(x) nchar(x) > 0, .)
dadaR2 = data.frame(title=dada2R1titles,denoiseR2=dada2denoisedR2)

dada2_denoise = left_join(dadaR1,dadaR2)

# fastq files
dada2FastqR1 = list.files(fastqFolder,pattern="R1_001.fastq.gz",full.names = TRUE)
dada2FastqR1titles = basename(dada2FastqR1) %>% strsplit(.,"_filt_S1_L001_R1_001.fastq.gz") %>% unlist(.) %>% strsplit(.,"_") %>% as.data.frame(.)
dada2FastqR1titles = dada2FastqR1titles[1,] %>% as.character(.)
dada2R1fastq = data.frame(title=dada2FastqR1titles,fastqR1=dada2FastqR1)

dada2FastqR2 = list.files(fastqFolder,pattern="R2_001.fastq.gz",full.names = TRUE)
dada2FastqR2titles = basename(dada2FastqR2) %>% strsplit(.,"_filt_S1_L001_R2_001.fastq.gz") %>% unlist(.) %>% strsplit(.,"_") %>% as.data.frame(.)
dada2FastqR2titles = dada2FastqR2titles[1,] %>% as.character(.)
dada2R2fastq = data.frame(title=dada2FastqR2titles,fastqR2=dada2FastqR2)

dada2_fastq = left_join(dada2R1fastq,dada2R2fastq)

dada2_files = left_join(dada2_fastq,dada2_denoise)

rm(dada2R2fastq,dada2FastqR2titles,dada2FastqR2,dada2R1fastq,dada2FastqR1titles,dada2FastqR1,dada2_denoise,dada2_fastq,dadaR2,dada2R2titles,dada2denoisedR2,dadaR1,dada2R1titles,dada2denoisedR1)

process_row <- function(row) {
  # Load the R1 data
  load(row['denoiseR1'])
  dadaFs <- dadaDenoised
  
  # Load the R2 data
  load(row['denoiseR2'])
  dadaRs <- dadaDenoised
  
  # Perform the merging operation
  print(paste("Now merging denoised R1 and R2 files from sample:", as.character(row['title'])))
  mergers = mergePairs(dadaFs, row['fastqR1'], dadaRs, row['fastqR2'], verbose = TRUE,minOverlap=5)
  
  return(mergers)
}

process_forward <- function(row) {
  # Load the R1 data
  load(row['denoiseR1'])
  dadaFs <- dadaDenoised
  
  # Perform the merging operation
  print(paste("Processing forward reads of:", as.character(row['title'])))
  forwardOnly = makeSequenceTable(dadaFs)
  
  return(forwardOnly)
}

# Apply the function to each row of the dataframe
mergers = apply(dada2_files, 1, process_row)

print("Sample merging completed.")
seqtab = makeSequenceTable(mergers) %>% as.data.frame(.)
seqtab2 = data.frame(SampleID=dada2_files$title,seqtab)

forwardOnly = apply(dada2_files, 1, process_forward)
seqtab3 = makeSequenceTable(forwardOnly) %>% as.data.frame(.)
seqtab4 = data.frame(SampleID=dada2_files$title,seqtab3)


write_tsv(seqtab2,file=paste0(theFolder,"/seqtab.tsv"),col_names=TRUE)
write_tsv(seqtab4,file=paste0(theFolder,"/seqtab_forwardOnly.tsv"),col_names=TRUE)
