#!/usr/bin/env Rscript
# compile_benchmarks.R
#
# Reads all Snakemake benchmark files from data/benchmarks/ and produces
# a compiled summary TSV for the manuscript.
#
# Snakemake benchmark columns:
#   s          — wall clock time in seconds
#   h:m:s      — wall clock time in h:m:s format
#   max_rss    — maximum resident set size (MB)
#   max_vms    — maximum virtual memory size (MB)
#   max_uss    — maximum unique set size (MB)
#   max_pss    — maximum proportional set size (MB)
#   io_in      — I/O read (MB)
#   io_out     — I/O written (MB)
#   mean_load  — mean CPU load (%)
#   cpu_time   — total CPU time (seconds)
#
# Usage:
#   Rscript --vanilla compile_benchmarks.R <output_file>

suppressMessages(library(dplyr))
suppressMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: compile_benchmarks.R <output_file>")
}
output_file <- args[1]

benchmark_dir <- "data/benchmarks"

# Find all benchmark files
files <- list.files(benchmark_dir, pattern = "\\.txt$", full.names = TRUE)

if (length(files) == 0) {
  stop("No benchmark files found in ", benchmark_dir)
}

message("Found ", length(files), " benchmark files.")

# Read and combine all benchmark files
benchmarks <- lapply(files, function(f) {
  df <- tryCatch(
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      message("Warning: could not read ", f, " — ", e$message)
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Extract rule name and wildcards from filename
  # Format: rule_name-wildcard1-wildcard2.txt
  fname <- tools::file_path_sans_ext(basename(f))

  df$file <- fname
  df
})

benchmarks <- bind_rows(Filter(Negate(is.null), benchmarks))

if (nrow(benchmarks) == 0) {
  stop("No valid benchmark data found.")
}

# Parse the filename into rule and wildcards
# Convention: first token before the first '-' that matches a known rule,
# or everything up to the first wildcard separator.
# Since rule names use underscores and wildcards are separated by '-',
# we split on the first '-' that follows the rule name pattern.
benchmarks <- benchmarks %>%
  mutate(
    # Extract rule: everything before the first wildcard
    # Rule names contain underscores; wildcards are joined with '-'
    rule = sub("^(.*?)-(batch|region|sample|[A-Z]|[a-z0-9]+[-].*)\\.?.*$", "\\1", file),
    # For files without wildcards (e.g., paired_taxonomy.txt), rule = file
    rule = ifelse(rule == file, file, rule),
    wildcards = ifelse(rule == file, "", sub(paste0("^", rule, "-?"), "", file))
  )

# Rename columns for clarity
if ("s" %in% colnames(benchmarks)) {
  benchmarks <- benchmarks %>%
    rename(
      wallclock_seconds = s,
      peak_memory_MB = max_rss,
      virtual_memory_MB = max_vms,
      cpu_seconds = cpu_time
    )
}

# Safe sd that returns NA instead of erroring when n < 2
sd_safe <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

# Create per-rule summary (aggregated across wildcards)
rule_summary <- benchmarks %>%
  group_by(pipeline_rule = rule) %>%
  summarise(
    n_files_processed       = n(),
    total_wallclock_seconds  = round(sum(wallclock_seconds, na.rm = TRUE), 1),
    mean_wallclock_seconds   = round(mean(wallclock_seconds, na.rm = TRUE), 1),
    sd_wallclock_seconds     = round(sd_safe(wallclock_seconds), 2),
    median_wallclock_seconds = round(median(wallclock_seconds, na.rm = TRUE), 1),
    min_wallclock_seconds    = round(min(wallclock_seconds, na.rm = TRUE), 1),
    max_wallclock_seconds    = round(max(wallclock_seconds, na.rm = TRUE), 1),
    mean_peak_memory_MB      = round(mean(peak_memory_MB, na.rm = TRUE), 1),
    sd_peak_memory_MB        = round(sd_safe(peak_memory_MB), 1),
    max_peak_memory_MB       = round(max(peak_memory_MB, na.rm = TRUE), 1),
    total_cpu_seconds        = round(sum(cpu_seconds, na.rm = TRUE), 1),
    mean_cpu_seconds         = round(mean(cpu_seconds, na.rm = TRUE), 1),
    sd_cpu_seconds           = round(sd_safe(cpu_seconds), 2),
    mean_cpu_load_pct        = round(mean(mean_load, na.rm = TRUE), 1),
    total_disk_read_MB       = round(sum(io_in, na.rm = TRUE), 1),
    total_disk_write_MB      = round(sum(io_out, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(total_wallclock_seconds))

# Write outputs
# Detailed per-job benchmarks
detail_file <- sub("\\.tsv$", "_detail.tsv", output_file)
benchmarks_out <- benchmarks %>%
  rename(
    pipeline_rule = rule,
    wallclock_hms = `h:m:s`,
    cpu_load_pct = mean_load,
    disk_read_MB = io_in,
    disk_write_MB = io_out
  ) %>%
  select(pipeline_rule, wildcards, wallclock_seconds, wallclock_hms,
         peak_memory_MB, virtual_memory_MB,
         cpu_seconds, cpu_load_pct, disk_read_MB, disk_write_MB) %>%
  arrange(pipeline_rule, wildcards)
write_tsv(benchmarks_out, detail_file)
message("Detailed benchmarks: ", detail_file)

# Per-rule summary
write_tsv(rule_summary, output_file)
message("Rule summary: ", output_file)

# Print summary to log
message("\n=== Benchmark Summary (by rule, sorted by total wall time) ===")
total_wall <- sum(rule_summary$total_wallclock_seconds)
total_cpu <- sum(rule_summary$total_cpu_seconds)
peak_rss <- max(rule_summary$max_peak_memory_MB)

# Helper: format seconds as human-readable duration
fmt_dur <- function(sec) {
  if (is.na(sec)) return("NA")
  if (sec < 60) return(sprintf("%.1fs", sec))
  if (sec < 3600) return(sprintf("%.1fm", sec / 60))
  sprintf("%.1fh", sec / 3600)
}

# Helper: format sd (show "n/a" when single job)
fmt_sd <- function(sd_val) {
  if (is.na(sd_val)) return("n/a")
  sprintf("%.2f", sd_val)
}

for (i in seq_len(nrow(rule_summary))) {
  r <- rule_summary[i, ]
  message(sprintf(
    "%-35s  files: %3d | wall: %s total, %s mean +/- %s sd [%s - %s] | RAM: %.0f MB mean, %.0f MB peak | CPU: %s total",
    r$pipeline_rule, r$n_files_processed,
    fmt_dur(r$total_wallclock_seconds), fmt_dur(r$mean_wallclock_seconds), fmt_sd(r$sd_wallclock_seconds),
    fmt_dur(r$min_wallclock_seconds), fmt_dur(r$max_wallclock_seconds),
    r$mean_peak_memory_MB, r$max_peak_memory_MB,
    fmt_dur(r$total_cpu_seconds)
  ))
}
message(sprintf("\nPipeline totals: wall=%s (%s)  CPU=%s (%s)  peak RSS=%.0f MB",
                fmt_dur(total_wall), sprintf("%.1f min", total_wall / 60),
                fmt_dur(total_cpu), sprintf("%.1f min", total_cpu / 60),
                peak_rss))
message("==============================================================")
