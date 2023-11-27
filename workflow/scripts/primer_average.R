suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))


option_list = list(
  make_option(c("-i", "--input"),  type="character", default=NULL, help="Input primer 1 OTU table", metavar="Input primer 1 OTU table"),
  make_option(c("-p", "--primer"), type="character", default=NULL, help="Input primer 2 OTU table", metavar="Input primer 2 OTU table"),
  make_option(c("-t", "--taxonomy"), type="character", default=NULL, help="Input taxonomy 1 table that has the OTU and the taxonomy", metavar="Input taxonomy 1 table that has the OTU and the taxonomy"),
  make_option(c("-t", "--taxa"),     type="character", default=NULL, help="Input taxonomy 2 table that has the OTU and the taxonomy", metavar="Input taxonomy 2 table that has the OTU and the taxonomy"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file", metavar="Output file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$database)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

primer1 = read_tsv(opt$input) %>% as.data.frame(.)
primer2 = read_tsv(opt$primer) %>% as.data.frame(.)
taxa1   = read_tsv(opt$taxonomy) %>% as.data.frame(.)
taxa2   = read_tsv(opt$taxa) %>% as.data.frame(.)

colnames(primer1)[1] = "OTU_ID"
colnames(primer2)[1] = "OTU_ID"

total_taxa = c(taxa1$Taxonomy,taxa2$Taxonomy) %>% unique(.)

test = matrix(nrow=length(total_taxa),ncol=dim(primer1)[2]) %>% as.data.frame(.)
colnames(test) = colnames(primer1)
test[,1] = total_taxa

primer1 = dplyr::left_join(primer1,taxa1) %>% .[,-which(colnames(.) == "Sequences")]
primer2 = dplyr::left_join(primer2,taxa2) %>% .[,-which(colnames(.) == "Sequences")]

heightOfTable = dim(test)[1]

for (x in 1:heightOfTable){ # Iterate over the rows
  currentTaxa = test[x,1]
  
  rowV13 = primer1[ which(primer1$Taxonomy == currentTaxa), ]
  rowV13$Taxonomy = NULL
  rowV13$OTU_ID = NULL
  rowV13 = colSums(rowV13) %>% t(.) %>% as.data.frame(.)
  
  rowV45 = primer2[ which(primer2$Taxonomy == currentTaxa), ]
  rowV45$Taxonomy = NULL
  rowV45$OTU_ID = NULL
  rowV45 = colSums(rowV45) %>% t(.) %>% as.data.frame(.)
  
  if (dim(rowV45)[1] == 0){
    test[x,] = rowV13 # If the rowV45 is empty, then add the V13 row
  } else if (dim(rowV13)[1] == 0){
    test[x,] = rowV45 # If the rowV13 is empty, then add the V45 row
  } else { # Meaning, the taxa exists in both tables
    
    V13count = rowV13 %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    V45count = rowV45 %>% as.numeric(.) %>% dplyr::na_if(0) %>% as.data.frame(.)
    
    meansOfRow = rowMeans(cbind(V13count,V45count),na.rm=TRUE) %>% round(.) # Meaning, if any of them has NA, then just report the normal number, then round everything
    meansOfRow[is.nan(meansOfRow)] = 0 # If na is present in both, then NaN will be reported, in this case, it means both rows had zeros, so just return it back to zero
    
    test[x,] = c(currentTaxa,meansOfRow)
    
  }
}

write_tsv(test,opt$output,col_names = TRUE)