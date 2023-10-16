suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="R1 object from denoise procedure", metavar="R1 object from denoise procedure"),
  make_option(c("-s", "--second"), type="character", default=NULL, help="R2 object from denoise procedure", metavar="R2 object from denoise procedure")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

R1denoised = opt$input
R2denoised = opt$second

load(R1denoised)

dadaFs = dadaDenoised
filtFs = theFiles

rm(theFiles, R1denoised, dadaDenoised, theFiles, DADA2errors)

load(R2denoised)

dadaRs = dadaDenoised
filtRs = theFiles

rm(theFiles, R1denoised, dadaDenoised, theFiles, DADA2errors)

print("Samples merging to begin…")
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
print("Sample merging completed.")
seqtab = makeSequenceTable(mergers)
seqtab = as.data.frame(cbind(rownames(seqtab),seqtab))

write_tsv(seqtab,file=paste0(dirname(inputFolder),"/seqtab.tsv"))

