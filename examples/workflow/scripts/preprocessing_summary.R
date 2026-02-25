#!/usr/bin/env Rscript
# preprocessing_summary.R
#
# Generates a summary table of the FAVABEAN preprocessing pipeline:
#   - Reads before/after quality filtering (per sample, per batch-region)
#   - Reads after chimera removal (per sample, per region)
#   - Number of ASVs per region
#   - Species-level assignment rate per region
#
# Usage:
#   Rscript --vanilla preprocessing_summary.R <output_file> [db_version]
#
# Arguments:
#   output_file  - Path for the per-region summary TSV
#   db_version   - (Optional) Versioned database name, e.g. "eHOMD_V15.22"
#                  Passed by Snakemake from the taxonomy_database config URL.
#                  Falls back to parsing the database name from ASV filenames.
#
# Expects the standard FAVABEAN directory layout under data/favabean/
# and reads files_info_Batches.csv for sample metadata.

suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: preprocessing_summary.R <output_file> [db_version]")
}
output_file <- args[1]
db_version  <- if (length(args) >= 2) args[2] else NULL

# ---- 1. Read sample metadata ----
mapping <- read.csv("data/files_info_Batches.csv", stringsAsFactors = FALSE)

# Build the expected filtered filename for each sample (matches dada2_1filterAndTrim.R)
# Input FASTQs are named {sample}_S{SampleNum}_L001_R1_001.fastq.gz after cutadapt.
# filterAndTrim splits on "_R1_" to get sample.names, then appends "_filt_S1_L001_R1_001.fastq.gz".
mapping$filt_name <- paste0(
  mapping$sample, "_S", mapping$SampleNum,
  "_L001_filt_S1_L001_R1_001.fastq.gz"
)

batch_regions <- unique(paste0(mapping$Batch_ID, "-", mapping$region))
regions <- unique(mapping$region)

# ---- 2. Collect filter statistics (reads before/after quality filtering) ----
filter_stats_list <- list()
for (br in batch_regions) {
  stats_file <- file.path("data/favabean", br, "dada2", "filter_stats.tsv")
  if (file.exists(stats_file)) {
    df <- read.delim(stats_file, stringsAsFactors = FALSE)
    df$batch_region <- br
    filter_stats_list[[br]] <- df
  } else {
    message("Warning: filter stats not found for ", br, " (", stats_file, ")")
  }
}

if (length(filter_stats_list) > 0) {
  filter_stats <- bind_rows(filter_stats_list)
  # The 'sample' column contains the filtered filename; join to mapping
  filter_stats <- filter_stats %>%
    left_join(
      mapping %>% select(filt_name, alias, region, Batch_ID),
      by = c("sample" = "filt_name")
    ) %>%
    # A sample sequenced across batches may appear in multiple filter_stats files.
    # Aggregate to one row per (alias, region) so downstream joins are 1:1.
    group_by(alias, region) %>%
    summarise(
      reads.in  = sum(reads.in,  na.rm = TRUE),
      reads.out = sum(reads.out, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  message("Warning: no filter stats files found. Pre-filter read counts will be NA.")
  # One row per unique (alias, region) — not per mapping row
  filter_stats <- mapping %>%
    distinct(alias, region) %>%
    mutate(reads.in = NA_real_, reads.out = NA_real_)
}

# ---- 3. Collect post-chimera read counts per region ----
chimera_stats_list <- list()
for (reg in regions) {
  robj_file <- file.path("data/favabean", paste0(reg, "_chimeraRemoved.RObjects"))
  if (file.exists(robj_file)) {
    load(robj_file)  # loads seqtab_noChim
    # rowSums = total reads per sample after chimera removal
    # ncol = number of unique ASVs after chimera removal
    per_sample <- data.frame(
      filt_name = rownames(seqtab_noChim),
      reads_after_dada2_5_remove_chimeras = rowSums(seqtab_noChim),
      stringsAsFactors = FALSE
    )
    per_sample$region <- reg
    per_sample$unique_ASVs_after_dada2_5_remove_chimeras <- ncol(seqtab_noChim)
    chimera_stats_list[[reg]] <- per_sample
    rm(seqtab_noChim)
  } else {
    message("Warning: chimera RObject not found for region ", reg)
  }
}

if (length(chimera_stats_list) > 0) {
  chimera_stats <- bind_rows(chimera_stats_list)
  # Join to mapping to get alias
  chimera_stats <- chimera_stats %>%
    left_join(
      mapping %>% select(filt_name, alias),
      by = "filt_name"
    ) %>%
    filter(!is.na(alias)) %>%
    # A sample sequenced across batches appears as separate rows in the merged
    # seqtab (mergeSequenceTables rbinds). Aggregate to one row per (alias, region)
    # to match what dada2_6condense.R does (group_by(alias) %>% summarize(sum)).
    group_by(alias, region) %>%
    summarise(
      reads_after_dada2_5_remove_chimeras  = sum(reads_after_dada2_5_remove_chimeras),
      unique_ASVs_after_dada2_5_remove_chimeras = first(unique_ASVs_after_dada2_5_remove_chimeras),
      .groups = "drop"
    )
} else {
  chimera_stats <- data.frame(
    alias = character(), region = character(),
    reads_after_dada2_5_remove_chimeras = numeric(),
    unique_ASVs_after_dada2_5_remove_chimeras = numeric()
  )
}

# ---- 4. Collect post-condense ASV counts per region ----
condense_stats_list <- list()
for (reg in regions) {
  condense_file <- file.path("data/favabean", paste0(reg, "_condense.tsv"))
  if (file.exists(condense_file)) {
    # Format: row 1 = "SampleIDs" + sample names, rows 2..N = sequence + counts
    condense <- read.delim(condense_file, header = FALSE, stringsAsFactors = FALSE)
    sample_names <- as.character(condense[1, -1])
    n_asvs_condensed <- nrow(condense) - 1
    # Per-sample total reads after condensing
    counts_mat <- condense[-1, -1, drop = FALSE]
    counts_mat <- apply(counts_mat, 2, as.numeric)
    if (!is.matrix(counts_mat)) counts_mat <- matrix(counts_mat, nrow = 1)
    per_sample <- data.frame(
      alias = sample_names,
      reads_after_dada2_6_condense = colSums(counts_mat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    per_sample$region <- reg
    per_sample$unique_ASVs_after_dada2_6_condense <- n_asvs_condensed
    condense_stats_list[[reg]] <- per_sample
  } else {
    message("Warning: condense file not found for region ", reg)
  }
}

if (length(condense_stats_list) > 0) {
  condense_stats <- bind_rows(condense_stats_list)
} else {
  condense_stats <- data.frame(
    alias = character(), region = character(),
    reads_after_dada2_6_condense = numeric(),
    unique_ASVs_after_dada2_6_condense = numeric()
  )
}

# ---- 5. Species-level assignment rate per region ----
taxonomy_stats_list <- list()
# Find all ASV files (they contain taxonomy column): {region}_{db}_ASV.tsv
asv_files <- list.files("data/favabean", pattern = "_ASV\\.tsv$", full.names = TRUE)
message(sprintf("Found %d ASV file(s) for taxonomy stats: %s",
                length(asv_files), paste(basename(asv_files), collapse = ", ")))
for (f in asv_files) {
  df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  # Extract region and database from filename: {region}_{db}_ASV.tsv
  fname <- basename(f)
  parts <- strsplit(sub("_ASV\\.tsv$", "", fname), "_")[[1]]
  reg <- parts[1]
  # Use versioned database name from config if available, otherwise parse filename
  db <- if (!is.null(db_version)) db_version else paste(parts[-1], collapse = "_")

  if ("taxonomy" %in% colnames(df)) {
    n_total <- nrow(df)
    # Species assigned: check both prefix-style (s__Species) and raw-style (7th rank)
    has_species <- vapply(df$taxonomy, function(tx) {
      # Prefix-style databases (SILVA, Greengenes, DADA2 output)
      if (grepl("s__[^;]+\\S", tx)) return(TRUE)
      # Raw-style databases (eHOMD): split on ";" and check 7th rank
      ranks <- trimws(unlist(strsplit(tx, ";")))
      ranks <- ranks[nchar(ranks) > 0]
      if (length(ranks) >= 7 && nchar(ranks[7]) > 0 && ranks[7] != ranks[6]) return(TRUE)
      FALSE
    }, logical(1), USE.NAMES = FALSE)
    n_species <- sum(has_species)
    species_rate <- round(n_species / n_total * 100, 1)
    taxonomy_stats_list[[fname]] <- data.frame(
      region = reg, taxonomy_database = db,
      unique_ASVs_after_dada2_7_assign_taxonomy = n_total,
      ASVs_assigned_to_species = n_species,
      pct_ASVs_assigned_to_species = species_rate,
      stringsAsFactors = FALSE
    )
  }
}

if (length(taxonomy_stats_list) > 0) {
  taxonomy_stats <- bind_rows(taxonomy_stats_list)
} else {
  taxonomy_stats <- data.frame(
    region = character(), taxonomy_database = character(),
    unique_ASVs_after_dada2_7_assign_taxonomy = numeric(),
    ASVs_assigned_to_species = numeric(),
    pct_ASVs_assigned_to_species = numeric()
  )
}

# ---- 5b. SIDLE reconstruction stats (when multiple regions exist) ----
#
# Two output tables are produced by SIDLE Rule 8:
#   - reconstructed_ASV.tsv         : accession-level (one row per ref clade)
#   - reconstructed_ASV_species.tsv : species-condensed (identical taxa summed)
# We report stats for both so the user can see the effect of condensation.
sidle_stats <- NULL
sidle_accession_file <- file.path("data/favabean", "reconstructed_ASV.tsv")
sidle_species_file   <- file.path("data/favabean", "reconstructed_ASV_species.tsv")
sidle_recon_tax      <- file.path("data/favabean/sidle", "reconstructed_taxonomy.tsv")

# Helper: detect whether a taxonomy string resolves to species
detect_species <- function(tx) {
  if (grepl("s__[^;]+\\S", tx)) return(TRUE)
  ranks <- trimws(unlist(strsplit(tx, ";")))
  ranks <- ranks[nchar(ranks) > 0]
  if (length(ranks) >= 7 && nchar(ranks[7]) > 0 && ranks[7] != ranks[6]) return(TRUE)
  FALSE
}

if (length(regions) > 1 && file.exists(sidle_accession_file)) {
  message("SIDLE reconstruction detected — collecting multi-region stats")

  # --- Accession-level stats ---
  acc_df <- read.delim(sidle_accession_file, stringsAsFactors = FALSE,
                       check.names = FALSE, comment.char = "")
  n_accession_features <- nrow(acc_df)
  acc_sample_cols <- setdiff(colnames(acc_df), c(colnames(acc_df)[1], "taxonomy"))
  recon_per_sample <- data.frame(
    alias = acc_sample_cols,
    reads_after_sidle_reconstruct = colSums(acc_df[, acc_sample_cols, drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Accession-level species rate
  acc_species_rate <- NA_real_
  if ("taxonomy" %in% colnames(acc_df)) {
    acc_has_sp <- vapply(acc_df$taxonomy, detect_species, logical(1), USE.NAMES = FALSE)
    acc_species_rate <- round(sum(acc_has_sp) / length(acc_has_sp) * 100, 1)
  }

  # --- Species-condensed stats ---
  n_species_features <- NA_integer_
  sp_species_rate    <- NA_real_
  if (file.exists(sidle_species_file)) {
    sp_df <- read.delim(sidle_species_file, stringsAsFactors = FALSE,
                        check.names = FALSE, comment.char = "")
    n_species_features <- nrow(sp_df)
    if ("taxonomy" %in% colnames(sp_df)) {
      sp_has_sp <- vapply(sp_df$taxonomy, detect_species, logical(1), USE.NAMES = FALSE)
      sp_species_rate <- round(sum(sp_has_sp) / length(sp_has_sp) * 100, 1)
    }
  }

  sidle_stats <- list(
    per_sample           = recon_per_sample,
    n_accession_features = n_accession_features,
    n_species_features   = n_species_features,
    acc_species_rate     = acc_species_rate,
    sp_species_rate      = sp_species_rate,
    n_regions            = length(regions)
  )
}

# ---- 6. Build per-sample summary ----
# All three tables now have at most one row per (alias, region),
# so these joins are guaranteed 1:1.
sample_summary <- filter_stats %>%
  left_join(
    chimera_stats %>% select(alias, region, reads_after_dada2_5_remove_chimeras),
    by = c("alias", "region")
  ) %>%
  left_join(
    condense_stats %>% select(alias, region, reads_after_dada2_6_condense),
    by = c("alias", "region")
  ) %>%
  rename(
    total_reads_raw = reads.in,
    reads_after_dada2_1_filter_and_trim = reads.out
  ) %>%
  arrange(region, alias)

# If SIDLE ran, add per-sample reconstructed counts
if (!is.null(sidle_stats)) {
  sample_summary <- sample_summary %>%
    left_join(sidle_stats$per_sample, by = "alias")
}

# ---- 7. Build per-region summary ----
# Helper: return NA if all values are NA, otherwise sum non-NA values
sum_or_na <- function(x) ifelse(all(is.na(x)), NA_real_, sum(x, na.rm = TRUE))

region_summary <- sample_summary %>%
  group_by(region) %>%
  summarise(
    n_samples = n(),
    total_reads_raw = sum_or_na(total_reads_raw),
    total_reads_after_dada2_1_filter_and_trim = sum_or_na(reads_after_dada2_1_filter_and_trim),
    total_reads_after_dada2_5_remove_chimeras = sum_or_na(reads_after_dada2_5_remove_chimeras),
    total_reads_after_dada2_6_condense = sum_or_na(reads_after_dada2_6_condense),
    .groups = "drop"
  ) %>%
  mutate(
    pct_retained_after_dada2_1_filter_and_trim = ifelse(
      !is.na(total_reads_raw) & total_reads_raw > 0,
      round(total_reads_after_dada2_1_filter_and_trim / total_reads_raw * 100, 1), NA_real_),
    pct_retained_after_dada2_5_remove_chimeras = ifelse(
      !is.na(total_reads_after_dada2_1_filter_and_trim) & total_reads_after_dada2_1_filter_and_trim > 0,
      round(total_reads_after_dada2_5_remove_chimeras / total_reads_after_dada2_1_filter_and_trim * 100, 1), NA_real_)
  ) %>%
  left_join(
    chimera_stats %>%
      select(region, unique_ASVs_after_dada2_5_remove_chimeras) %>% distinct(),
    by = "region"
  ) %>%
  left_join(
    condense_stats %>%
      select(region, unique_ASVs_after_dada2_6_condense) %>% distinct(),
    by = "region"
  ) %>%
  left_join(taxonomy_stats, by = "region")

# ---- 8. Write outputs ----
# Per-sample detail
sample_output <- sub("\\.tsv$", "_per_sample.tsv", output_file)
write_tsv(sample_summary, sample_output)
message("Per-sample summary: ", sample_output)

# Per-region summary
write_tsv(region_summary, output_file)
message("Per-region summary: ", output_file)

# Print to log as well
message("\n=== Preprocessing Summary (per region) ===")
for (i in seq_len(nrow(region_summary))) {
  r <- region_summary[i, ]
  message(sprintf(
    "Region: %s | Samples: %d | Raw reads: %s | After dada2_1_filter_and_trim: %s (%s%%) | After dada2_5_remove_chimeras: %s (%s%%) | ASVs after chimera: %d (after condense: %d) | Species assigned: %.1f%%",
    r$region, r$n_samples,
    format(r$total_reads_raw, big.mark = ","),
    format(r$total_reads_after_dada2_1_filter_and_trim, big.mark = ","),
    ifelse(is.na(r$pct_retained_after_dada2_1_filter_and_trim), "NA",
           sprintf("%.1f", r$pct_retained_after_dada2_1_filter_and_trim)),
    format(r$total_reads_after_dada2_5_remove_chimeras, big.mark = ","),
    ifelse(is.na(r$pct_retained_after_dada2_5_remove_chimeras), "NA",
           sprintf("%.1f", r$pct_retained_after_dada2_5_remove_chimeras)),
    r$unique_ASVs_after_dada2_5_remove_chimeras, r$unique_ASVs_after_dada2_6_condense,
    ifelse(is.na(r$pct_ASVs_assigned_to_species), 0, r$pct_ASVs_assigned_to_species)
  ))
}

# SIDLE multi-region reconstruction summary
if (!is.null(sidle_stats)) {
  message("\n--- SIDLE Multi-Region Reconstruction ---")
  message(sprintf("  Regions combined: %d (%s)",
                  sidle_stats$n_regions, paste(regions, collapse = ", ")))
  message(sprintf("  Unique taxa after sidle_reconstruct: %d | Species assigned: %s%%",
                  sidle_stats$n_accession_features,
                  ifelse(is.na(sidle_stats$acc_species_rate), "NA",
                         sprintf("%.1f", sidle_stats$acc_species_rate))))
  message(sprintf("  Unique taxa after sidle_to_biom: %s | Species assigned: %s%%",
                  ifelse(is.na(sidle_stats$n_species_features), "NA",
                         as.character(sidle_stats$n_species_features)),
                  ifelse(is.na(sidle_stats$sp_species_rate), "NA",
                         sprintf("%.1f", sidle_stats$sp_species_rate))))

  # Write SIDLE-specific summary file alongside the main summary
  sidle_summary <- data.frame(
    taxonomy_database                   = if (!is.null(db_version)) db_version else "unknown",
    n_regions_combined                  = sidle_stats$n_regions,
    regions                             = paste(regions, collapse = ","),
    unique_taxa_after_sidle_reconstruct             = sidle_stats$n_accession_features,
    pct_species_after_sidle_reconstruct      = sidle_stats$acc_species_rate,
    unique_taxa_after_sidle_to_biom        = sidle_stats$n_species_features,
    pct_species_after_sidle_to_biom = sidle_stats$sp_species_rate,
    stringsAsFactors = FALSE
  )
  sidle_output <- sub("\\.tsv$", "_sidle.tsv", output_file)
  write_tsv(sidle_summary, sidle_output)
  message("SIDLE summary: ", sidle_output)
}

message("==========================================")
