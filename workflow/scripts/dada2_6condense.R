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
rownames(seqtab_noChim3) = seqtab_noChim3$alias
seqtab_noChim3$alias = NULL

print("Starting collapse of ASVs…")

collapseNoMismatchOptimized <- function(seqtab, minOverlap = 20, orderBy = "abundance", 
                                        identicalOnly = FALSE, vec = TRUE, band = -1, 
                                        verbose = FALSE) {
  if (verbose) message("Starting optimized collapseNoMismatch function.")
  
  # Remove duplicate sequences
  dupes <- duplicated(colnames(seqtab))
  if (any(dupes)) {
    st <- seqtab[, !dupes, drop = FALSE]
    for (i in which(dupes)) {
      sq <- colnames(seqtab)[[i]]
      st[, sq] <- st[, sq] + seqtab[, i]
    }
    seqtab <- st
  }
  if (identicalOnly) {
    return(seqtab)
  }
  
  # Get unique sequences and sort by abundance
  unqs.srt <- sort(colSums(seqtab), decreasing = TRUE)
  seqs <- names(unqs.srt)
  
  # Precompute prefixes
  prefixes <- substr(seqs, 1, minOverlap)
  names(prefixes) <- seqs  # Map sequences to their prefixes
  
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
  if (verbose) message("Output ", ncol(collapsed), " collapsed sequences out of ", 
                       ncol(seqtab), " input sequences.")
  
  return(collapsed)
}




optimize_minOverlap <- function(seqtab, min_overlap_range) {
  max_collapsed <- 0
  optimal_minOverlap <- min_overlap_range[1]
  
  for (minOverlap in min_overlap_range) {
    print(paste("Testing overlap area size:",minOverlap))
    result <- collapseNoMismatchOptimized(
      seqtab,
      minOverlap = minOverlap,
      orderBy = "abundance",
      identicalOnly = FALSE,
      vec = TRUE,
      band = -1,
      verbose = FALSE
    )
    num_collapsed <- ncol(result)
    
    if (num_collapsed < max_collapsed) {
      # Stop if the number of collapsed sequences decreases
      break
    } else {
      max_collapsed <- num_collapsed
      optimal_minOverlap <- minOverlap
    }
  }
  
  return(optimal_minOverlap)
}

# Define the range of minOverlap values
min_overlap_range <- seq(5, 100, by = 5)

# Find the optimal minOverlap
optimal_minOverlap <- optimize_minOverlap(as.matrix(seqtab_noChim3), min_overlap_range)

# Ensure minOverlap_values is a numeric vector
minOverlap_values <- as.numeric(minOverlap_values)

# Use the Kneedle algorithm to find the elbow point
elbow_point <- inflection::findiplist(
  x = minOverlap_values,
  y = num_collapsed,
  index = 1
)


# Extract the optimal index using matrix indexing
optimal_index <- elbow_point["ESE", "j1"]
optimal_index <- as.integer(optimal_index)

# Check if the index is valid
if (is.na(optimal_index) || optimal_index < 1 || optimal_index > length(minOverlap_values)) {
  stop("Optimal index is invalid.")
}

# Get the optimal minOverlap value
optimal_minOverlap <- minOverlap_values[optimal_index]



# Number of bootstrap samples
n_boot <- 1000

# Initialize vectors to store inflection points
inflection_points_ESE <- numeric(n_boot)
inflection_points_EDE <- numeric(n_boot)

set.seed(123) # For reproducibility

for (i in 1:n_boot) {
  # Resample indices with replacement
  resample_indices <- sample(length(minOverlap_values), replace = TRUE)
  boot_minOverlap <- minOverlap_values[resample_indices]
  boot_num_collapsed <- num_collapsed[resample_indices]
  
  # Ensure data is sorted for the inflection function
  sorted_indices <- order(boot_minOverlap)
  boot_minOverlap <- boot_minOverlap[sorted_indices]
  boot_num_collapsed <- boot_num_collapsed[sorted_indices]
  
  # Apply inflection point detection
  boot_elbow_point <- inflection::findiplist(
    x = boot_minOverlap,
    y = boot_num_collapsed,
    index = 1
  )
  
  # Extract optimal indices
  optimal_index_ESE <- as.integer(boot_elbow_point["ESE", "j1"])
  optimal_index_EDE <- as.integer(boot_elbow_point["EDE", "j1"])
  
  # Handle possible NA values
  if (!is.na(optimal_index_ESE) && optimal_index_ESE >= 1 && optimal_index_ESE <= length(boot_minOverlap)) {
    inflection_points_ESE[i] <- boot_minOverlap[optimal_index_ESE]
  } else {
    inflection_points_ESE[i] <- NA
  }
  
  if (!is.na(optimal_index_EDE) && optimal_index_EDE >= 1 && optimal_index_EDE <= length(boot_minOverlap)) {
    inflection_points_EDE[i] <- boot_minOverlap[optimal_index_EDE]
  } else {
    inflection_points_EDE[i] <- NA
  }
}

# Remove NA values
inflection_points_ESE <- inflection_points_ESE[!is.na(inflection_points_ESE)]
inflection_points_EDE <- inflection_points_EDE[!is.na(inflection_points_EDE)]

# Calculate confidence intervals
ci_ESE <- quantile(inflection_points_ESE, probs = c(0.025, 0.975))
ci_EDE <- quantile(inflection_points_EDE, probs = c(0.025, 0.975))

# Function to check if two intervals overlap
intervals_overlap <- function(ci1, ci2) {
  return(!(ci1[2] < ci2[1] || ci2[2] < ci1[1]))
}

# Check if confidence intervals overlap
overlap <- intervals_overlap(ci_ESE, ci_EDE)

# Calculate widths of confidence intervals
width_ESE <- ci_ESE[2] - ci_ESE[1]
width_EDE <- ci_EDE[2] - ci_EDE[1]

cat("Confidence Interval Widths - ESE:", width_ESE, ", EDE:", width_EDE, "\n")

# Initialize selected estimator
selected_estimator <- NULL

if (overlap) {
  # Confidence intervals overlap
  cat("Confidence intervals overlap.\n")
  # Choose the estimator with the smaller confidence interval
  if (width_ESE < width_EDE) {
    selected_estimator <- "ESE"
  } else {
    selected_estimator <- "EDE"
  }
} else {
  # Confidence intervals do not overlap
  cat("Confidence intervals do not overlap.\n")
  # Choose the estimator with the confidence interval entirely within the other
  if (ci_ESE[1] > ci_EDE[1] && ci_ESE[2] < ci_EDE[2]) {
    selected_estimator <- "ESE"
  } else if (ci_EDE[1] > ci_ESE[1] && ci_EDE[2] < ci_ESE[2]) {
    selected_estimator <- "EDE"
  } else {
    # If neither is entirely within the other, choose the one with the narrower confidence interval
    if (width_ESE < width_EDE) {
      selected_estimator <- "ESE"
    } else {
      selected_estimator <- "EDE"
    }
  }
}

cat("Exact Slope Estimator (ESE). Exact Derivative Estimator (EDE).")
cat("ESE: Estimates the inflection point by analyzing changes in the slope of the data.")
cat("EDE: Estimates the inflection point by examining changes in the curvature of the data")
cat("Selected Estimator:", selected_estimator, "\n")


# Based on the selected estimator, compute the optimal minOverlap
if (selected_estimator == "ESE") {
  # Use the median of the ESE bootstrap estimates
  optimal_minOverlap <- median(inflection_points_ESE)
} else if (selected_estimator == "EDE") {
  # Use the median of the EDE bootstrap estimates
  optimal_minOverlap <- median(inflection_points_EDE)
} else {
  stop("Estimator selection failed.")
}

cat("Optimal minOverlap determined by", selected_estimator, "is:", optimal_minOverlap, "\n")


# Ensure optimal_minOverlap is one of the values in minOverlap_values
optimal_minOverlap <- minOverlap_values[which.min(abs(minOverlap_values - optimal_minOverlap))]

# Run the function with the selected minOverlap
condensed_table <- collapseNoMismatchOptimized(
  as.matrix(seqtab_noChim3),
  minOverlap = optimal_minOverlap,
  orderBy = "abundance",
  identicalOnly = FALSE,
  vec = TRUE,
  band = -1,
  verbose = TRUE)


SampleIDs = rownames(condensed_table)
condensed_table = as.data.frame(condensed_table)
condensed_table = cbind(SampleIDs,condensed_table)
condensed_table = t(condensed_table) %>% as.data.frame(.)
condensed_table = cbind(rownames(condensed_table),condensed_table)

write_tsv(condensed_table,opt$output,col_names = FALSE)