suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library("dada2"))

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

if (any(is.na(seqtab_noChim3$alias))==TRUE){seqtab_noChim3 = seqtab_noChim3[-which(is.na(seqtab_noChim3$alias)),]}
rownames(seqtab_noChim3) = seqtab_noChim3$alias
seqtab_noChim3$alias = NULL

collapseNoMismatchOptimized <- function(seqtab, minOverlap = 20, orderBy = "abundance", 
                                        identicalOnly = FALSE, vec = TRUE, band = -1, 
                                        verbose = FALSE) {
  if (verbose) message("Starting optimized collapseNoMismatch function.")
  
  # Get unique sequences and sort by abundance
  unqs.srt <- sort(colSums(seqtab), decreasing = TRUE)
  seqs <- names(unqs.srt)
  
  # Precompute prefixes
  prefixes <- substr(seqs, 1, minOverlap)
  names(prefixes) <- seqs  # Map sequences to their prefixes
  message("Prefixes computed. Starting grouping…")
  # Group sequences by their prefixes
  prefix_groups <- split(seqs, prefixes)
  
  # Create an environment to store sequences and their counts
  seq_env <- new.env(hash = TRUE, parent = emptyenv())
  
  if (verbose) message("Processing sequences...")
  
  for (group_prefix in names(prefix_groups)) {
    group_seqs <- prefix_groups[[group_prefix]]
    group_seqtab <- seqtab[, group_seqs, drop = FALSE]
    
    # Sort sequences in the group by abundance
    group_abundance <- colSums(group_seqtab)
    group_seqs_sorted <- names(sort(group_abundance, decreasing = TRUE))
    
    rep_seq <- group_seqs_sorted[1]
    rep_counts <- group_seqtab[, rep_seq]
    
    for (seq in group_seqs_sorted[-1]) {
      if (nwhamming(seq, rep_seq, vec = vec, band = band) == 0) {
        # Sequences are identical after alignment; collapse counts
        rep_counts <- rep_counts + group_seqtab[, seq]
      } else {
        # Add new sequence to the environment
        assign(seq, group_seqtab[, seq], envir = seq_env)
      }
    }
    
    # Store the representative sequence and its counts
    assign(rep_seq, rep_counts, envir = seq_env)
  }
  
  # Collect sequences and counts from the environment
  collapsed_seqs <- ls(seq_env)
  collapsed_counts <- mget(collapsed_seqs, envir = seq_env)
  
  # Convert counts to a matrix
  collapsed <- do.call(cbind, collapsed_counts)
  colnames(collapsed) <- collapsed_seqs
  rownames(collapsed) <- rownames(seqtab)
  
  # Order sequences if needed
  if (!is.null(orderBy)) {
    if (orderBy == "abundance") {
      ord <- order(colSums(collapsed), decreasing = TRUE)
      collapsed <- collapsed[, ord, drop = FALSE]
    } else if (orderBy == "nsamples") {
      ord <- order(colSums(collapsed > 0), decreasing = TRUE)
      collapsed <- collapsed[, ord, drop = FALSE]
    }
  }
  
  if (verbose) message("Finished processing sequences.")
  if (verbose) message("Output: ",(ncol(seqtab)-ncol(collapsed)), " sequences were collapsed into other sequences. Result is ", ncol(collapsed), " collapsed sequences out of ", 
                       ncol(seqtab), " input sequences.")
  
  return(collapsed)
}

# Find the optimal minOverlap and collect minOverlap_values and num_collapsed
print("Starting collapse of ASVs…")
optimization_results <- collapseNoMismatchOptimized(as.matrix(seqtab_noChim3), verbose=TRUE)

SampleIDs = rownames(optimization_results)
condensed_table = as.data.frame(optimization_results)
condensed_table = cbind(SampleIDs,condensed_table)
condensed_table = t(condensed_table) %>% as.data.frame(.)
condensed_table = cbind(rownames(condensed_table),condensed_table)

write_tsv(condensed_table,opt$output,col_names = FALSE)
