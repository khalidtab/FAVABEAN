#!/usr/bin/env Rscript
# seqtab_to_fasta.R
# Convert DADA2 condensed output into SIDLE-compatible inputs:
#   1. FASTA file with stable ASV IDs (ASV_1, ASV_2, ...)
#   2. Count TSV with ASV IDs as row names, samples as columns
#
# Input:  {region}_condense.tsv (transposed seqtab from dada2_6condense.R)
#         Format: no header, col 1 = row IDs. Row 1 = "SampleIDs" + sample names.
#         Rows 2..N = DNA sequence + abundance counts per sample.
#
# Output: {region}_rep_seqs.fasta   — representative sequences
#         {region}_asv_counts.tsv   — count table with ASV IDs
#         {region}_seqlengths.txt   — sequence length stats (min, Q1, median, Q3, max)
#                                     first line = recommended trim length (Q1)
#
# Usage:
#   Rscript seqtab_to_fasta.R <input_condense.tsv> <output_fasta> <output_counts_tsv> <output_seqlengths>

suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript seqtab_to_fasta.R <condense.tsv> <out.fasta> <out_counts.tsv> <out_seqlengths>")
}

input_file       <- args[1]
output_fasta     <- args[2]
output_counts    <- args[3]
output_seqlengths <- args[4]

# Read the transposed condense table (no header, no col names)
raw <- read_tsv(input_file, col_names = FALSE, show_col_types = FALSE)

# Row 1 contains "SampleIDs" followed by sample names
sample_names <- as.character(raw[1, -1])

# Rows 2..N: col 1 = DNA sequence, cols 2..M = counts
sequences <- as.character(raw[-1, 1, drop = TRUE])
counts_mat <- as.data.frame(raw[-1, -1])
colnames(counts_mat) <- sample_names

# Convert counts to numeric
counts_mat <- as.data.frame(lapply(counts_mat, as.numeric))

# Assign stable ASV IDs
asv_ids <- paste0("ASV_", seq_along(sequences))
rownames(counts_mat) <- asv_ids

# Write FASTA
fasta_lines <- character(length(sequences) * 2)
for (i in seq_along(sequences)) {
  fasta_lines[(i - 1) * 2 + 1] <- paste0(">", asv_ids[i])
  fasta_lines[(i - 1) * 2 + 2] <- sequences[i]
}
writeLines(fasta_lines, output_fasta)

# Write count table (ASV IDs as first column)
out_df <- cbind(data.frame(ASV_ID = asv_ids, stringsAsFactors = FALSE), counts_mat)
write_tsv(out_df, output_counts)

# Compute sequence length statistics
seq_lengths <- nchar(sequences)
len_stats <- quantile(seq_lengths, probs = c(0, 0.25, 0.5, 0.75, 1))
names(len_stats) <- c("min", "Q1", "median", "Q3", "max")

# Write trim length file:
#   Line 1: recommended trim length (Q1 — retains >= 75% of sequences)
#   Line 2+: full stats for logging
trim_length <- as.integer(len_stats["Q1"])
stats_lines <- c(
  trim_length,
  sprintf("min=%d Q1=%d median=%d Q3=%d max=%d n_asvs=%d",
          len_stats["min"], len_stats["Q1"], len_stats["median"],
          len_stats["Q3"], len_stats["max"], length(seq_lengths))
)
writeLines(as.character(stats_lines), output_seqlengths)

message(sprintf("Converted %d ASVs across %d samples", length(sequences), length(sample_names)))
message(sprintf("  FASTA:  %s", output_fasta))
message(sprintf("  Counts: %s", output_counts))
message(sprintf("  Sequence lengths: min=%d, Q1=%d, median=%d, Q3=%d, max=%d",
                len_stats["min"], len_stats["Q1"], len_stats["median"],
                len_stats["Q3"], len_stats["max"]))
message(sprintf("  Recommended SIDLE trim length: %d bp", trim_length))
