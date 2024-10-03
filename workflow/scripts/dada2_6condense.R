suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library("dada2"))
suppressMessages(library("inflection"))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file after chimera removal", metavar="input file after chimera removal"),
  make_option(c("-r", "--region"), type="character", default=NULL, help="input region", metavar="input region"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="output file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



inputFile = opt$input
cores = opt$cores

load(inputFile)

seqtab_noChim = as.data.frame(seqtab_noChim)

mapping = read.csv("data/files_info_Batches.csv")
mapping$theFiles = paste0(mapping$sample,"_S",mapping$SampleNum,"_L001_filt_S1_L001_R1_001.fastq.gz")
mapping2 = data.frame(thefiles = mapping$theFiles, alias = mapping$alias)

seqtab_noChim2 = cbind(rownames(seqtab_noChim),seqtab_noChim)
colnames(seqtab_noChim2)[1] = "thefiles"

seqtab_noChim3 = dplyr::left_join(seqtab_noChim2,mapping2)
seqtab_noChim3$thefiles = NULL
print("Samples sequenced multiple times (if any) will now be merged…")
seqtab_noChim3 = seqtab_noChim3 %>% group_by(alias) %>% summarize(across(everything(), sum)) %>% as.data.frame(.)
print("Sample combining (if any) done.")
rownames(seqtab_noChim3) = seqtab_noChim3$alias
seqtab_noChim3$alias = NULL

print("Starting collapse of ASVs…")
condensed_table = collapseNoMismatch(as.matrix(seqtab_noChim3), minOverlap = 20, orderBy = "abundance", identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)

SampleIDs = rownames(condensed_table)
condensed_table = as.data.frame(condensed_table)
condensed_table = cbind(SampleIDs,condensed_table)
condensed_table = t(condensed_table) %>% as.data.frame(.)
condensed_table = cbind(rownames(condensed_table),condensed_table)

write_tsv(condensed_table,opt$output,col_names = FALSE)