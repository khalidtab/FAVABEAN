suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))
suppressMessages(library("tictoc"))

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

rm(theFiles, dadaDenoised)

load(R2denoised)

dadaRs = dadaDenoised
filtRs = theFiles

rm(theFiles, dadaDenoised)

print("Samples merging to begin…")
tic()
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
toc()
print("Sample merging completed.")
seqtab = makeSequenceTable(mergers)

write_tsv(as.data.frame(cbind(rownames(seqtab),seqtab)),file=paste0(dirname(inputFolder),"/seqtab.tsv"),col_names=TRUE)

