suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("jsonlite"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Figaro json file", metavar="JSON file from Figaro, which shows you the DADA2 parameter options"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="The number of cores to use in this script", metavar="The number of cores to use in this script"),
  make_option(c("-f", "--folder"), type="character", default=NULL, help="Folder with the fastq files to be read by DADA2", metavar="Folder of fastq"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Folder for the DADA2 trimmed and filtered fastq files", metavar="Folder for the DADA2 trimmed and filtered fastq files"),
  make_option(c("-p", "--option"), type="character", default=NULL, help="The parameter of interest to use in DADA2", metavar="DADA2 parameter")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

jsonInput = opt$input %>% jsonlite::fromJSON(.)
# Split the max expected error column into the two values
splitExpectedEE = as.data.frame(jsonInput$maxExpectedError) %>% as.matrix(.) %>% t(.)
colnames(splitExpectedEE) = c("R1EE","R2EE")
splitTrimPosition = as.data.frame(jsonInput$trimPosition) %>% as.matrix(.) %>% t(.)
colnames(splitTrimPosition) = c("R1Trim","R2Trim")
readRetentionPercent = jsonInput$readRetentionPercent
score = jsonInput$score

myJSONTable  =  cbind(splitExpectedEE,splitTrimPosition,readRetentionPercent,score) %>% as.data.frame(.)


inputFolder = opt$folder
outputFolder = opt$output
DADA2Options = opt$option 

cores = as.numeric(opt$cores)

if (DADA2Options == "highest_coverage"){

  most_coverage = myJSONTable %>% .[which(.$readRetentionPercent == max(.$readRetentionPercent)),]
  result_row = most_coverage[which.max(most_coverage$score),] # If there are still ties, find row with maximum score
  } 

if (DADA2Options == "least_errors"){

  min_R1EE_rows = myJSONTable %>% .[.$R1EE == min(.$R1EE),] # Find rows with minimum R1EE
  min_R2EE_rows <- min_R1EE_rows[min_R1EE_rows$R2EE == min(min_R1EE_rows$R2EE),] # If there are ties, find rows with minimum R2EE
  result_row = min_R2EE_rows[which.max(min_R2EE_rows$score),] # If there are still ties, find row with maximum score
  }

print("Parameters for DADA2 has been chosen.")

fnFs = sort(list.files(inputFolder, pattern="R1", full.names = TRUE))
fnRs = sort(list.files(inputFolder, pattern="R2", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_"), `[`, 1)

filtFs <- file.path(paste0(outputFolder), "filtered", paste0(sample.names, "_filt_S1_L001_R1_001.fastq.gz"))
filtRs <- file.path(paste0(outputFolder), "filtered", paste0(sample.names, "_filt_S1_L001_R2_001.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print("Files detected, will perform filtering and trimming on the fastq files.")

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    truncLen=c(result_row$R1Trim,result_row$R2Trim),
                    maxN=0,
                    maxEE=c(result_row$R1EE,result_row$R2EE), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=cores, matchIDs=TRUE)

