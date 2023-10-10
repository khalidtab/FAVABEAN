suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("jsonlite"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Figaro json file", metavar="JSON file from Figaro, which shows you the DADA2 parameter options"),
  make_option(c("-f", "--folder"), type="character", default=NULL, help="Folder with the fastq files to be read by DADA2", metavar="Folder of fastq"),
  make_option(c("-p", "--option"), type="character", default=NULL, help="The parameter of interest to use in DADA2", metavar="DADA2 parameter")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

jsonInput = opt$input %>% jsonlite::fromJSON(.) # jsonInput = "/Users/khaled/Desktop/v13_V45/cutadapt/run5-V13/figaro/trimParameters.json" %>% jsonlite::fromJSON(.)
# Split the max expected error column into the two values
splitExpectedEE = as.data.frame(jsonInput$maxExpectedError) %>% as.matrix(.) %>% t(.)
colnames(splitExpectedEE) = c("R1EE","R2EE")
splitTrimPosition = as.data.frame(jsonInput$trimPosition) %>% as.matrix(.) %>% t(.)
colnames(splitTrimPosition) = c("R1Trim","R2Trim")
readRetentionPercent = jsonInput$readRetentionPercent
score = jsonInput$score

myJSONTable  =  cbind(splitExpectedEE,splitTrimPosition,readRetentionPercent,score) %>% as.data.frame(.)


inputFolder = opt$folder # inputFolder = "/Users/khaled/Desktop/v13_V45/cutadapt/run5-V13/cutadapt"
DADA2Options = opt$option # DADA2options = "most_coverage" # options: most_coverage, least_errors, best score

if (DADA2Options == "most_coverage"){
  most_coverage = myJSONTable %>% .[which(.$readRetentionPercent == max(.$readRetentionPercent)),]
  result_row = most_coverage[which.max(most_coverage$score),] # If there are still ties, find row with maximum score
  
} else if (DADA2Options == "least_errors"){

  min_R1EE_rows = myJSONTable %>% .[.$R1EE == min(.$R1EE),] # Find rows with minimum R1EE
  min_R2EE_rows <- min_R1EE_rows[min_R1EE_rows$R2EE == min(min_R1EE_rows$R2EE),] # If there are ties, find rows with minimum R2EE
  result_row = min_R2EE_rows[which.max(min_R2EE_rows$score),] # If there are still ties, find row with maximum score
  }

fnFs = sort(list.files(inputFolder, pattern="R1", full.names = TRUE))
fnRs = sort(list.files(inputFolder, pattern="R2", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_"), `[`, 1)
filtFs <- file.path(paste0(inputFolder,"/dada2"), "filtered", paste0(sample.names, "_filt_S1_L001_R1_001.fastq.gz"))
filtRs <- file.path(paste0(inputFolder,"/dada2"), "filtered", paste0(sample.names, "_filt_S1_L001_R2_001.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    truncLen=c(result_row$R1Trim,result_row$R2Trim),
                    maxN=0,
                    maxEE=c(result_row$R1EE,result_row$R2EE), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE, matchIDs=TRUE)

errF = learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
errR = learnErrors(filtRs, multithread=TRUE, randomize=TRUE)

# Dereplicated
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)

# Denoising the sequences
dadaFs = dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs = dada(derepRs, err=errR, multithread=TRUE, pool=TRUE) # This resulted in barely any sequences, so we are only going to use

mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab = makeSequenceTable(mergers)
seqtab = as.data.frame(cbind(rownames(seqtab),seqtab))

write_tsv(seqtab,file=paste0(inputFolder,"/seqtab.tsv"))

