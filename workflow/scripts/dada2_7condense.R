suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input folder with all the seqtab files", metavar="input folder with all the seqtab files")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



inputFolder = opt$input

seqtabs = list.files(inputFolder, pattern="seqtab.tsv", full.names = TRUE, recursive = TRUE)

combinedTables = data.frame(V1=as.character(),
                      sequence=as.character(),
                      abundance=as.numeric())

for (x in seqtabs){
currentTable =  read_tsv(x) %>% as.data.frame(.) 
currentTable = currentTable %>% pivot_longer(.,cols=2:dim(currentTable)[2], names_to = "sequence",values_to = "abundance") %>% as.data.frame(.)
combinedTables = rbind(combinedTables,currentTable)
}



seqtab_combined = pivot_wider(combinedTables, names_from = "sequence", values_from = "abundance")
seqtab_combined = as.data.frame(seqtab_combined)
row.names(seqtab_combined) = seqtab_combined[,1]
seqtab_combined[,1] = NULL
seqtab_combined = as.matrix(seqtab_combined)

print("All tables combined")

seqtab_noChim =  removeBimeraDenovo(seqtab_combined, method = "pooled", multithread =  cores, verbose=TRUE)

print("Chimeras removed")

condensed_table = collapseNoMismatch(seqtab_noChim, minOverlap = 20, orderBy = "abundance", identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)
SampleIDs = rownames(condensed_table)
condensed_table = as.data.frame(condensed_table)
condensed_table = cbind(SampleIDs,condensed_table)

region = strsplit(inputFolder,"/dada2/") %>% as.character(.) %>% basename(.) %>% strsplit(.,"-")
region = region[[1]][2]

pathForCondensed_table = strsplit(inputFolder,"favabean")[[1]][1]



write_tsv(condensed_table,
          paste0(pathForCondensed_table,"/favabean/",region,"_condensed_ASVs.tsv"))
