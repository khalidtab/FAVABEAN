suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))


option_list = list(
  make_option(c("-d", "--database"), type="character", default=NULL, help="Which database to use", metavar="Which database to use")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$database)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}



myDatabase = opt$database


myFiles = list.files("data/favabean",pattern=paste0(myDatabase,"_taxonomy.tsv"),full.names = TRUE)

myFile1 = read_tsv(myFiles[1])
myFile2 = read_tsv(myFiles[2])

total_taxa = c(myFile1$taxonomy,myFile2$taxonomy) %>% unique(.)

test = matrix(nrow=length(total_taxa),ncol=dim(myFile2)[2]) %>% as.data.frame(.)
colnames(test) = colnames(myFile2)
test[,1] = total_taxa

heightOfTable = dim(test)[1]

for (x in 1:heightOfTable){ # Iterate over the rows
  currentTaxa = test[x,1]
  
  rowV13 = myFile1[ which(myFile1$taxonomy == currentTaxa), ]
  rowV13 = rowV13[,-1]
  rowV13$taxonomy = NULL
  rowV13 = colSums(rowV13) %>% t(.) %>% as.data.frame(.)
  
  rowV45 = myFile2[ which(myFile2$taxonomy == currentTaxa), ]
  rowV45 = rowV45[,-1]
  rowV45$taxonomy = NULL
  rowV45 = colSums(rowV45) %>% t(.) %>% as.data.frame(.)
  
  if (dim(rowV45)[1] == 0){
    test[x,] = rowV13 # If the rowV45 is empty, then add the V13 row
  } else if (dim(rowV13)[1] == 0){
    test[x,] = rowV45 # If the rowV13 is empty, then add the V45 row
  } else { # Meaning, the taxa exists in both tables
    
    V13count = rowV13 %>% as.character(.) %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    V45count = rowV45 %>% as.character(.) %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    
    meansOfRow = rowMeans(cbind(V13count,V45count),na.rm=TRUE) %>% round(.) # Meaning, if any of them has NA, then just report the normal number, then round everything
    meansOfRow[is.nan(meansOfRow)] = 0 # If na is present in both, then NaN will be reported, in this case, it means both rows had zeros, so just return it back to zero
    
    test[x,] = c(currentTaxa,meansOfRow,currentTaxa)
    
  }
}

write_tsv(test,paste0("data/favabean/primerAveraged_",myDatabase,".tsv"),col_names = TRUE)