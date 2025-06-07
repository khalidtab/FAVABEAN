suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))

theFiles = list.files("data/favabean/","_taxonomy.tsv",full.names = TRUE)

readAndCondense = function(primer_path,taxa_path){
  primer = read_tsv(primer_path,show_col_types = FALSE) %>% as.data.frame(.)
  taxa   = read_tsv(taxa_path,show_col_types = FALSE) %>% as.data.frame(.)  
  colnames(primer)[1] = "OTU_ID"
  primer_joined = dplyr::left_join(primer,taxa) %>% .[,-which(colnames(.) %in% c("Sequences"))]
  
  message("Condensing the sequences based on the taxonomic rank.")
  primer_joined_condensed = primer_joined %>% select(-OTU_ID) %>% group_by(Taxonomy) %>% summarize(across(everything(), sum)) %>% as.data.frame(.)

  return(list(primer_joined   =  as.data.frame(primer_joined),
              primer_condensed=as.data.frame(primer_joined_condensed)))
  }

opt = NULL
opt$output = "data/favabean/primer_averaged.tsv"


if (length(theFiles) == 1){message("There is only one primer. Skipping primer averaging and finalizing biom file creation")}else{
  
  
message("First primer.")
primer1 = read_tsv(theFiles[1]) %>% as.data.frame(.)


message("Second primer.")
primer2 = read_tsv(theFiles[2]) %>% as.data.frame(.)

sampleNames = unique(c(colnames(primer1), colnames(primer2)))
sampleNames = sampleNames[-which(sampleNames %in% c("taxonomy","#SampleID"))]

sharedNonSharedSamples = (c(colnames(primer1),colnames(primer2)))
shared = sharedNonSharedSamples[duplicated(sharedNonSharedSamples)]
shared = shared[-which(shared %in% c("taxonomy","#SampleID"))]

notShared = sampleNames[-which(shared %in% sampleNames)]

if(length(notShared) > 0){
  warning("The following samples were found in only one primer. Printing information on the samples, and will continue primer averaging without them.")
  print(notShared)
  print(paste0("Primer1: ", basename(theFiles[1]),
        "    primer2: ",basename(theFiles[2])))
}

primer1 = primer1 %>% select("#SampleID",shared)
primer2 = primer2 %>% select("#SampleID",shared)

taxonomy = unique(c(primer1[,1],primer2[,1]))

test = matrix(nrow=length(taxonomy),ncol=1+length(shared)) %>% as.data.frame(.)
colnames(test) = c("OTU_ID",shared)
test[,1] = taxonomy
heightOfTable = dim(test)[1]

for (x in 1:heightOfTable){ # Iterate over the rows
  currentTaxa = test[x,1]
  
  rowV13 = primer1[ which(primer1[,1] == currentTaxa), ]
  rowV13$taxonomy = NULL
  rowV13$`#SampleID` = NULL
  
  rowV45 = primer2[ which(primer2[,1] == currentTaxa), ]
  rowV45$taxonomy = NULL
  rowV45$`#SampleID` = NULL
  
  if (is.na(sum(as.numeric(rowV45)))){
    test[x,c(2,3)] = rowV13 # If the rowV45 is empty, then add the V13 row
  } else if (is.na(sum(as.numeric(rowV13)))){
    test[x,c(2,3)] = rowV45 # If the rowV13 is empty, then add the V45 row
  } else { # Meaning, the taxa exists in both tables
    
    V13count = rowV13 %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    V45count = rowV45 %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    
    meansOfRow = rowMeans(cbind(V13count,V45count),na.rm=TRUE) %>% round(.) # Meaning, if any of them has NA, then just report the normal number, then round everything
    meansOfRow[is.nan(meansOfRow)] = 0 # If na is present in both, then NaN will be reported, in this case, it means both rows had zeros, so just return it back to zero
    
    test[x,] = c(currentTaxa,meansOfRow)
    
  }
}

test[,1] = taxonomy
test[, -1] <- lapply(test[, -1], function(x) as.numeric(as.character(x)))

test[is.na(test)] = 0

write_tsv(test,"data/favabean/primer_averaged.tsv",col_names = TRUE)
}


