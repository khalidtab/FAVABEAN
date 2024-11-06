suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library("dada2"))


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file after ASV condensation", metavar="input file after ASV condensation"),
  make_option(c("-d", "--database"), type="character", default=NULL, help="input database to use for taxonomy assignemnt", metavar="input database"),
  make_option(c("-s", "--species"), type="character", default=NULL, help="input database for species assignment to use for taxonomy assignemnt", metavar="input species database"),
  make_option(c("-c", "--cores"), type="character", default=NULL, help="Number of cores to use", metavar="Number of cores to use"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="output file"),
  make_option(c("-t", "--taxonomy"), type="character", default=NULL, help="output taxonomy file", metavar="output taxonomy file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

inputFile = opt$input
myDatabase = opt$database
mySpecies = opt$species
myCores = as.numeric(opt$cores)
myOutput = opt$output
myCondensed = opt$taxonomy

inputMatrix = read_tsv(inputFile) %>% as.data.frame(.)
ASVs = inputMatrix$SampleIDs


print("Beginning assignment up to genus level.")
taxa3 = dada2::assignTaxonomy(ASVs, myDatabase,verbose=TRUE, tryRC = TRUE, multithread = TRUE)

print("Taxonomy assignment up to genus level is completed.")

print("Beginning assignment up to species level.")
taxa4 = dada2::assignSpecies(taxa3,mySpecies,tryRC = TRUE)

print("Taxonomy assignment up to species level is completed.")

taxa3 = cbind(rownames(taxa3),taxa3) %>% as.data.frame(.)
taxa4 = cbind(rownames(taxa4),taxa4) %>% as.data.frame(.)

taxa5 = dplyr::left_join(taxa3,taxa4)
taxa5[is.na(taxa5)] = ""
taxonomy = paste0("k__",taxa5$Kingdom,"; p__",taxa5$Phylum,"; c__",taxa5$Class,"; o__",taxa5$Order,"; f__",taxa5$Family,"; g__",taxa5$Genus,"; s__",taxa5$Species)


inputMatrix2 = cbind(inputMatrix,taxonomy)
OTU_IDs = rownames(inputMatrix2)
inputMatrix3 = cbind(paste0("OTU",OTU_IDs),inputMatrix2)

colnames(inputMatrix3)[1] = "#SampleID"
taxonomy = data.frame(OTU_ID = inputMatrix3$`#SampleID`, Sequences = inputMatrix3$SampleIDs, Taxonomy = inputMatrix3$taxonomy)

inputMatrix3 = inputMatrix3 %>% .[,-which(colnames(.) %in% c("SampleIDs","taxonomy"))]

# Now let's create the different iterations of the files needed
readAndCondense = function(primer_path,taxa_path){
  primer = as.data.frame(primer_path)
  taxa   = as.data.frame(taxa_path)
  colnames(primer)[1] = "OTU_ID"
  primer_joined = dplyr::left_join(primer,taxa) %>% .[,-which(colnames(.) %in% c("Sequences"))]
  
  colnames(primer_joined)[length(colnames(primer_joined))] = "taxonomy"
  
  message("Condensing the sequences based on the taxonomic rank.")
  primer_joined_condensed = primer_joined %>% select(-OTU_ID) %>% group_by(taxonomy) %>% summarize(across(everything(), sum)) %>% as.data.frame(.)
  
  
  return(list(primer_joined   =  as.data.frame(primer_joined),
              primer_condensed=as.data.frame(primer_joined_condensed)))
}

myIterations = readAndCondense(inputMatrix3,taxonomy)
inputMatrix4 = left_join(inputMatrix3,taxonomy,by = join_by(`#SampleID` == OTU_ID)) %>% select(-Sequences)
colnames(inputMatrix4)[length(colnames(inputMatrix4))] = "taxonomy"
write_tsv(inputMatrix4,myOutput)

colnames(myIterations$primer_condensed)[1] = "#SampleID"

write_tsv(myIterations[["primer_condensed"]],myCondensed)

