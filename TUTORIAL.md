# FAVABEAN Tutorial: From Raw FASTQs to Taxonomy Tables

This tutorial walks you through a complete run of the FAVABEAN pipeline using the bundled example dataset. By the end, you will have ASV count tables with taxonomy assignments for each amplicon region, a multi-region reconstructed table combining all regions via the SMURF EM algorithm (SIDLE), and pipeline summary reports — all generated automatically from raw paired-end FASTQ files.

The example dataset included in `examples/` contains 2 biological samples (P1A1, P1A2) sequenced across 2 amplicon regions (27F targeting V1–V3, and V3V5 targeting V3–V5) on 2 Illumina sequencing runs. This results in 16 FASTQ files (2 samples × 2 regions × 2 runs × 2 read directions), each containing 2,000 reads.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [1. Setting Up Your Data](#1-setting-up-your-data)
  - [1.1 Prepare files_info.csv](#11-prepare-files_infocsv)
  - [1.2 Configure favabean.yaml](#12-configure-favabenyaml)
- [2. Running the Pipeline](#2-running-the-pipeline)
  - [2.1 Quick Start (Automated)](#21-quick-start-automated)
  - [2.2 Step-by-Step (Manual)](#22-step-by-step-manual)
- [3. Understanding the Pipeline Steps](#3-understanding-the-pipeline-steps)
- [4. Output Files — What You Get](#4-output-files--what-you-get)
  - [4.1 Primary Outputs (in favabean/)](#41-primary-outputs-in-favabean)
  - [4.2 Per-Region ASV and Taxonomy Tables](#42-per-region-asv-and-taxonomy-tables)
  - [4.3 Summary Reports](#43-summary-reports)
  - [4.4 Intermediate Files (Batch Folders)](#44-intermediate-files-batch-folders)
- [5. Interpreting the Summary Reports](#5-interpreting-the-summary-reports)
  - [5.1 Preprocessing Summary (Per-Region)](#51-preprocessing-summary-per-region)
  - [5.2 Preprocessing Summary (Per-Sample)](#52-preprocessing-summary-per-sample)
  - [5.3 Benchmark Summary](#53-benchmark-summary)
- [6. Advanced: SIDLE Multi-Region Reconstruction](#6-advanced-sidle-multi-region-reconstruction)
  - [6.1 When Does SIDLE Run?](#61-when-does-sidle-run)
  - [6.2 SIDLE Configuration](#62-sidle-configuration)
  - [6.3 SIDLE Outputs](#63-sidle-outputs)
  - [6.4 SIDLE Summary Report](#64-sidle-summary-report)
  - [6.5 SIDLE Intermediate Files](#65-sidle-intermediate-files)
- [7. Using Outputs with FALAPhyl](#7-using-outputs-with-falaphyl)
- [8. Troubleshooting](#8-troubleshooting)
- [9. FAQ](#9-faq)

---

## Prerequisites

You need one of the following setups:

**Option A — Docker (recommended):**

```bash
docker pull khalidtab/favabean:latest
```

**Option B — Local installation:**

Snakemake (≥ 7.0) and Conda (or Mamba) installed on your system. The pipeline creates its own Conda environments for each tool automatically.


## 1. Setting Up Your Data

### 1.1 Prepare files_info.csv

The `files_info.csv` file tells FAVABEAN about your samples. Each row maps a FASTQ pair to its sample name, primer sequences, amplicon region, and expected amplicon length. Here is the example dataset's `files_info.csv`:

```
sample,fastq1,fastq2,primer_5,primer_3,region,Batch_ID,expected length
P1A1,P1A1_27F_S1346_L001_R1_001.fastq.gz,P1A1_27F_S1346_L001_R2_001.fastq.gz,ACAC...CTCAG,GTGA...GCTG,27F,,500
P1A1,P1A1_27F_S3065_L001_R1_001.fastq.gz,P1A1_27F_S3065_L001_R2_001.fastq.gz,ACAC...CTCAG,GTGA...GCTG,27F,,500
P1A1,P1A1_V3V5_S1755_L001_R1_001.fastq.gz,...,V3V5,,500
...
```

The columns are:

| Column | Description |
|--------|-------------|
| `sample` | Biological sample name. Same sample sequenced multiple times should have the same name — FAVABEAN merges them automatically. |
| `fastq1` | Filename of the forward (R1) FASTQ file, relative to the `data/` directory. |
| `fastq2` | Filename of the reverse (R2) FASTQ file. |
| `primer_5` | Full forward primer sequence, including any Illumina adapter prefix. IUPAC degenerate bases (Y, R, W, etc.) are supported. |
| `primer_3` | Full reverse primer sequence, including any Illumina adapter prefix. |
| `region` | Amplicon region name (e.g., `27F`, `V3V5`, `V4`). Used to group samples for processing. |
| `Batch_ID` | Leave **empty** — FAVABEAN auto-detects sequencing batches from the FASTQ headers. |
| `expected length` | Expected amplicon insert length in base pairs (used by Figaro to optimize trim parameters). |


The way that FAVABEAN detects the Batch_ID is through reading the heading in the FASTQ files. Here is how it looks like in P1A1_27F_S1346_L001_R1_001:

> @VH01519:74:AAFHWY5M5:1:1101:35803:1000 1:N:0:GTCTGATACG+ACAGACTCTA

And here is how it looks like in P1A1_27F_S1346_L001_R1_001.fastq.gz:

> @VH01519:71:AAFHMVTM5:1:1101:37432:1000 1:N:0:GTCTGATACG+ACAGACTCTA

Since the `sample` is the same despite it being sequenced in two different sequencing runs, DADA2 errors should be learned for each separately. When you run the `start`, you will see instructions similar to this:

```
====================
Welcome to FAVABEAN!
====================


Copying the template files in your directly if they don't exist (environments.txt, files_info.csv, and favabean.yaml)'
If this is the first time you are running the pipeline on this dataset, you will need to initialize your 'files_info.csv' file.
Please write 'conda activate' and copy the biom.yaml environment link from FAVABEAN_environments.txt file

Example:
 conda activate <<biom environment such as : .snakemake/conda/007a7beaa3353a33c938a5a0e57be4ff_>> 

Then, execute the following command
 Rscript --vanilla workflow/scripts/batch.R


Example command(s):
snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete --scheduler greedy 
Please copy and paste the above command(s) or type your own command below.
```

The instructions, as they say, will require you to initialize the `files_info.csv` file. This basically means it will look at all your fastq files, and create the batches necessary. Once you do that, the pipeline will automatically learn the DADA2 errors for them separately then merge the sequences together after ASVs are generated. The `files_info_Batches.csv` will look like this

```
sample,fastq1,fastq2,primer_5,primer_3,region,Batch_ID,expected.length,alias,SampleNum
P1A1,P1A1_27F_S1346_L001_R1_001.fastq.gz,P1A1_27F_S1346_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAKRGTTYGATYNTGGCTCAG,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTACGTNTBACCGCDGCTGCTG,27F,Batch1,500,P1A1,1
P1A1_1,P1A1_27F_S3065_L001_R1_001.fastq.gz,P1A1_27F_S3065_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAKRGTTYGATYNTGGCTCAG,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTACGTNTBACCGCDGCTGCTG,27F,Batch2,500,P1A1,2
P1A1_2,P1A1_V3V5_S1755_L001_R1_001.fastq.gz,P1A1_V3V5_S1755_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGYCAGCMGCCGCGGTAA,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGYCAATTYMTTTRAGTTT,V3V5,Batch2,500,P1A1,3
P1A1_3,P1A1_V3V5_S36_L001_R1_001.fastq.gz,P1A1_V3V5_S36_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGYCAGCMGCCGCGGTAA,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGYCAATTYMTTTRAGTTT,V3V5,Batch1,500,P1A1,4
P1A2,P1A2_27F_S1353_L001_R1_001.fastq.gz,P1A2_27F_S1353_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAKRGTTYGATYNTGGCTCAG,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTACGTNTBACCGCDGCTGCTG,27F,Batch1,500,P1A2,5
P1A2_1,P1A2_27F_S3072_L001_R1_001.fastq.gz,P1A2_27F_S3072_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAKRGTTYGATYNTGGCTCAG,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTACGTNTBACCGCDGCTGCTG,27F,Batch2,500,P1A2,6
P1A2_2,P1A2_V3V5_S1761_L001_R1_001.fastq.gz,P1A2_V3V5_S1761_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGYCAGCMGCCGCGGTAA,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGYCAATTYMTTTRAGTTT,V3V5,Batch2,500,P1A2,7
P1A2_3,P1A2_V3V5_S42_L001_R1_001.fastq.gz,P1A2_V3V5_S42_L001_R2_001.fastq.gz,ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTGYCAGCMGCCGCGGTAA,GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGYCAATTYMTTTRAGTTT,V3V5,Batch1,500,P1A2,8
```

Note that the primer sequences should include the Illumina adapter prefix if it is present in your reads. FAVABEAN's SIDLE module automatically strips known Illumina adapter sequences when extracting the biological primer for in-silico PCR.

If you have only one amplicon region, the pipeline runs the standard DADA2 workflow. If you have two or more regions, SIDLE multi-region reconstruction is triggered automatically.

### 1.2 Configure favabean.yaml

The `favabean.yaml` file controls pipeline behavior. The example config (shown below) works well for most datasets:

```yaml
# Minimum sequence length after primer removal
initial_filter:
     - 200

# How to determine trim length for Figaro
# Options: default, max_len, Q1, Q2, Q3, or a specific number (e.g., 100)
trim_param:
     - default

# Which Figaro parameter set to use for DADA2
# Options: highest_coverage, lowest_errors
figaro:
    - highest_coverage

# SIDLE multi-region reconstruction parameters (only used with 2+ regions)
sidle:
  max_mismatch: 2
  primer_max_mismatch: 1
  min_counts: 1000
  min_abund: 1.0e-5
  region_normalize:
    - average

# Reference database for taxonomy assignment
taxonomy_database:
  eHOMD:
    use: True
    url: "https://zenodo.org/records/7380292/files/eHOMD_RefSeq_dada2_V15.22.fasta.gz"
    species: "https://zenodo.org/records/7380292/files/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz"
```

**Configuration options explained:**

`initial_filter` — Sequences shorter than this (in bp) after primer removal are discarded. 200 bp is a good default to remove primer dimers and very short fragments.

`trim_param` — Controls how the uniform trim length is chosen before running Figaro. The `default` option uses the 25th percentile (Q1) if it is within 10% of the median (Q2), otherwise uses Q2. You can also specify a fixed number.

`figaro` — Figaro evaluates many possible trim positions and scores them. `highest_coverage` picks the trim position that retains the most reads; `lowest_errors` picks the one with the fewest expected errors.

`taxonomy_database` — The pipeline downloads the reference database automatically on first run. You can switch to SILVA by setting `eHOMD.use: False` and `silva.use: True`:

```yaml
taxonomy_database:
  eHOMD:
    use: False
    url: "..."
    species: "..."
  silva:
    use: True
    url: "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
    species: "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz"
```


## 2. Running the Pipeline

### 2.1 Quick Start (Automated)

The fastest way to try FAVABEAN is to run the bundled example script. From the FAVABEAN project root:

```bash
bash examples/run_example.sh
```

This script copies the example FASTQ files, `files_info.csv`, and `favabean.yaml` into `data/`, then runs the full pipeline. On a modern laptop with 4 cores, the example dataset completes in roughly 5–10 minutes (longer on first run while Conda environments are created if you did not use the Docker environment).

### 2.2 Step-by-Step (Manual)

**With Docker:**

```bash
docker run --rm -it -v ~/path/to/your/data:/data khalidtab/favabean:latest start
```

This drops you into the container. Then:

```bash
snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 2 --rerun-incomplete --scheduler greedy
```

**Without Docker:**

Place your FASTQ files, `files_info.csv`, and `favabean.yaml` inside the `data/` directory of the FAVABEAN project, then from the project root:

```bash
snakemake paired_taxonomy \
    --snakefile workflow/snakefile \
    --use-conda \
    --cores all \
    --keep-going \
    --retries 2 \
    --rerun-incomplete \
    --scheduler greedy
```

The `--use-conda` flag is **required** — the pipeline creates separate Conda environments for each tool (cutadapt, DADA2, Figaro, SIDLE, etc.). On first run, Snakemake downloads and installs these environments, which can take 10–20 minutes. Subsequent runs reuse the cached environments.


## 3. Understanding the Pipeline Steps

FAVABEAN processes your data through the following stages. Each stage corresponds to a named Snakemake rule, and the rule name appears in the output summaries so you can trace exactly where each number comes from.

```
Raw FASTQs
    │
    ▼
trim_primers_with_cutadapt --- Remove primer and adapter sequences
    │
    ▼
sequence_length_stats -------- Calculate read length distribution
    │
    ▼
filter_length_for_figaro ----- Trim all reads to uniform length
    │
    ▼
estimate_trim_params_w_figaro - Find optimal DADA2 trim parameters
    │
    ▼
dada2_1_filter_and_trim ------ Quality filter and trim reads
    │
    ▼
dada2_2_learn_errors -------- Learn sequencing error rates (R1 and R2 separately)
    │
    ▼
dada2_3_denoise ------------- Denoise reads into ASVs (R1 and R2 separately)
    │
    ▼
dada2_4_merge_paired_ends ---- Merge R1 and R2 into full-length ASVs
    │
    ▼
dada2_5_remove_chimeras ------ Remove chimeric sequences (de novo)
    │
    ▼
dada2_6_condense ------------ Merge ASV tables across sequencing batches
    │
    ▼
dada2_7_assign_taxonomy ------ Assign taxonomy using reference database
    │
    ▼
Per-Region ASV + Taxonomy tables (.tsv and .biom)
    │
    ├-- (single region: you're done!)
    │
    └-- (2+ regions: SIDLE reconstruction runs automatically)
            │
            ▼
        SIDLE: Reconstruct a unified count table across all regions
            │
            ▼
        Reconstructed ASV tables (accession-level and species-condensed)
```

Steps that run **per batch-region combination** (e.g., Batch1-27F, Batch2-V3V5): 
- `trim_primers_with_cutadapt`
-  `sequence_length_stats`
- `filter_length_for_figaro`
- `estimate_trim_params_w_figaro`
- `dada2_1_filter_and_trim`
- `dada2_2_learn_errors`
- `dada2_3_denoise`
- `dada2_4_merge_paired_ends`.

Steps that run **per region** (combining all batches): 
- `dada2_5_remove_chimeras`
- `dada2_6_condense`
- `dada2_7_assign_taxonomy`

Steps that run **once** (combining all regions): SIDLE reconstruction, final summary generation.


## 4. Output Files — What You Get

After a successful run, your output directory (`data/favabean/` inside Docker, or the `favabean/` folder in the project root) contains the following files.

### 4.1 Primary Outputs (in favabean/)

These are the files you will use for downstream analysis:

| File | Description |
|------|-------------|
| `{region}_{db}_ASV.tsv` | Per-region ASV count table with taxonomy. Rows are ASVs, columns are samples. |
| `{region}_{db}_ASV.biom` | Same data in BIOM format (for tools like QIIME2, phyloseq, FALAPhyl). |
| `{region}_{db}_taxonomy.tsv` | Per-region taxonomy-only table. Rows are unique taxonomy strings, columns are summed counts. |
| `{region}_{db}_taxonomy.biom` | Same data in BIOM format. |
| `reconstructed_ASV.tsv` | Multi-region reconstructed count table (accession-level). One row per reference sequence. |
| `reconstructed_ASV.biom` | Same in BIOM format. |
| `reconstructed_ASV_species.tsv` | Species-condensed reconstruction. Rows with identical taxonomy are summed. |
| `reconstructed_ASV_species.biom` | Same in BIOM format. |
| `preprocessing_summary.tsv` | Per-region pipeline statistics. |
| `preprocessing_summary_per_sample.tsv` | Per-sample read counts at each pipeline stage. |
| `preprocessing_summary_sidle.tsv` | SIDLE-specific reconstruction statistics. |
| `benchmark_summary.tsv` | Computation time and memory usage per pipeline rule. |
| `benchmark_summary_detail.tsv` | Per-file detailed benchmarks. |

For the example dataset (2 regions: 27F, V3V5; database: eHOMD), the primary outputs are:

```
favabean/
├-- 27F_eHOMD_ASV.tsv
├-- 27F_eHOMD_ASV.biom
├-- 27F_eHOMD_taxonomy.tsv
├-- 27F_eHOMD_taxonomy.biom
├-- V3V5_eHOMD_ASV.tsv
├-- V3V5_eHOMD_ASV.biom
├-- V3V5_eHOMD_taxonomy.tsv
├-- V3V5_eHOMD_taxonomy.biom
├-- reconstructed_ASV.tsv
├-- reconstructed_ASV.biom
├-- reconstructed_ASV_species.tsv
├-- reconstructed_ASV_species.biom
├-- preprocessing_summary.tsv
├-- preprocessing_summary_per_sample.tsv
├-- preprocessing_summary_sidle.tsv
├-- benchmark_summary.tsv
└-- benchmark_summary_detail.tsv
```

### 4.2 Per-Region ASV and Taxonomy Tables

**ASV table** (`27F_eHOMD_ASV.tsv`):

```
#SampleID    P1A1    P1A2    taxonomy
ASV1         200     203     k__Bacteria; p__Firmicutes; ...; g__Streptococcus; s__
ASV2         125     144     k__Bacteria; p__Firmicutes; ...; g__Streptococcus; s__
ASV3         130     136     k__Bacteria; p__Firmicutes; ...; g__Streptococcus; s__
```

Each row is a unique ASV (amplicon sequence variant). The count columns show how many reads from each sample mapped to that ASV. The `taxonomy` column contains the taxonomic classification from kingdom through species.

**Taxonomy table** (`27F_eHOMD_taxonomy.tsv`):

```
#SampleID                                                              P1A1    P1A2
k__Bacteria; p__Actinobacteria; ...; g__Actinomyces; s__              248     186
k__Bacteria; p__Actinobacteria; ...; g__Actinomyces; s__massiliensis  15      0
```

This table groups all ASVs that share the same taxonomy string and sums their counts. This is useful when you want to work at the genus or species level rather than the individual ASV level.

### 4.3 Summary Reports

See [Section 5](#5-interpreting-the-summary-reports) for detailed interpretation.

### 4.4 Intermediate Files (Batch Folders)

Each batch-region combination gets its own folder (e.g., `favabean/Batch1-27F/`) containing intermediate files from each pipeline step:

```
favabean/Batch1-27F/
├-- cutadapt/                      # Primer-trimmed FASTQs
├-- initialFilt_discarded/         # Reads that were too short after trimming
├-- dada2/                         # Quality-filtered FASTQs
│   ├-- filter_stats.tsv           # Reads before/after quality filtering
│   └-- filtered/                  # Filtered FASTQ files
├-- figaro/
│   └-- trimParameters.json        # Optimal trim parameters chosen by Figaro
├-- dada2_DADA2Errors-R1.RData     # Learned error model (forward reads)
├-- dada2_DADA2Errors-R2.RData     # Learned error model (reverse reads)
├-- dada2_DADA2Denoise-R1.RData    # Denoised ASVs (forward)
├-- dada2_DADA2Denoise-R2.RData    # Denoised ASVs (reverse)
├-- seqtab.tsv                     # Merged paired-end ASV table
└-- stats_R1_lengths.txt           # Read length distribution
```

Additionally, per-region files appear directly in `favabean/`:

```
favabean/
├-- 27F_chimeraRemoved.RObjects    # Post-chimera ASV table (R object)
├-- 27F_condense.tsv               # Batch-merged ASV count matrix
├-- V3V5_chimeraRemoved.RObjects
└-- V3V5_condense.tsv
```

These intermediate files are useful for debugging but are not needed for downstream analysis.


## 5. Interpreting the Summary Reports

### 5.1 Preprocessing Summary (Per-Region)

`preprocessing_summary.tsv` gives you a high-level view of what happened to your reads at each pipeline step. Here is the example output:

| region | n_samples | total_reads_raw | total_reads_after\_dada2_1\_filter\_and\_trim | total_reads_after\_dada2_5\_remove\_chimeras | total_reads_after\_dada2_6\_condense | pct_retained_after\_dada2_1\_filter\_and\_trim | pct_retained_after\_dada2_5\_remove\_chimeras | unique_ASVs\_after\_dada2_5\_remove\_chimeras | unique_ASVs\_after\_dada2_6\_condense | taxonomy_database | unique_ASVs\_after\_dada2_7\_assign\_taxonomy | ASVs_assigned\_to\_species | pct_ASVs_assigned\_to\_species |
|--------|-----------|-----------------|------|------|------|------|------|-----|-----|------|-----|-----|------|
| 27F | 2 | 8000 | 7241 | 3693 | 3693 | 90.5 | 51.0 | 89 | 89 | eHOMD_V15.22 | 89 | 89 | 100 |
| V3V5 | 2 | 8000 | 7306 | 6374 | 6374 | 91.3 | 87.2 | 61 | 61 | eHOMD_V15.22 | 61 | 61 | 100 |

**Column names reference the pipeline rule that produced each value**, making it easy to trace where each number comes from:

- `total_reads_raw` — Total reads across all samples before any processing.
- `total_reads_after_dada2_1_filter_and_trim` — Reads surviving quality filtering.
- `pct_retained_after_dada2_1_filter_and_trim` — Percentage of raw reads that passed quality filtering. In the example, ~90% of reads passed. If this drops below 50%, consider relaxing your quality thresholds.
- `total_reads_after_dada2_5_remove_chimeras` — Reads remaining after chimera removal. The 27F region shows only 51% retention, which is normal for longer amplicons where chimera formation is more common.
- `unique_ASVs_after_dada2_5_remove_chimeras` — Number of unique ASV sequences after chimera filtering but before batch merging. If you sequenced on multiple runs, this number reflects the combined (but not deduplicated) count.
- `unique_ASVs_after_dada2_6_condense` — Number of unique ASVs after merging samples across batches. This can be the same as the chimera count (as in the example, where both batches produced the same ASVs) or lower if batches had overlapping ASVs that were merged.
- `unique_ASVs_after_dada2_7_assign_taxonomy` — Number of ASVs that received a taxonomy assignment. Should equal the condense count unless the reference database cannot classify some sequences.
- `pct_ASVs_assigned_to_species` — Percentage of ASVs resolved to species level. 100% in this example because eHOMD is a curated oral database with species annotations.

### 5.2 Preprocessing Summary (Per-Sample)

`preprocessing_summary_per_sample.tsv` tracks read counts for each individual sample through each step:

| alias | region | total_reads_raw | reads_after\_dada2_1\_filter\_and\_trim | reads_after\_dada2_5\_remove\_chimeras | reads_after\_dada2_6\_condense | reads_after\_sidle\_reconstruct |
|-------|--------|-----------------|------|------|------|------|
| P1A1 | 27F | 4000 | 3642 | 1975 | 1975 | 2084 |
| P1A2 | 27F | 4000 | 3599 | 1718 | 1718 | 2266 |
| P1A1 | V3V5 | 4000 | 3669 | 3016 | 3016 | 2084 |
| P1A2 | V3V5 | 4000 | 3637 | 3358 | 3358 | 2266 |

The `reads_after_sidle_reconstruct` column appears only when SIDLE ran (multiple regions). Note that these values are **per sample across all regions combined** — the SMURF EM algorithm reconstructs a single unified abundance table, so each sample gets one total count.

It is normal for the SIDLE reconstructed counts to differ from the per-region counts. The EM algorithm redistributes ambiguous reads using cross-region evidence, so counts can go up or down relative to any single region.

### 5.3 Benchmark Summary

`benchmark_summary.tsv` reports computation time and memory usage for each pipeline rule:

| pipeline_rule | n_files_processed | total_wallclock_seconds | mean_wallclock_seconds | mean_peak_memory_MB | max_peak_memory_MB |
|---|---|---|---|---|---|
| dada2_7_assign_taxonomy-27F-eHOMD | 1 | 19.9 | 19.9 | 767.0 | 767.0 |
| dada2_1_filter_and_trim | 2 | 16.2 | 8.1 | 724.1 | 725.5 |
| dada2_2_learn_errors-Batch2-27F | 2 | 14.0 | 7.0 | 749.3 | 778.7 |

This tells you which rules took the longest and used the most memory. For example, taxonomy assignment loaded the full eHOMD reference database into memory, using ~767 MB of RAM. The `n_files_processed` column tells you how many parallel jobs ran for that rule — useful for understanding the load on your system.

The `benchmark_summary_detail.tsv` file provides per-file granularity, showing the exact timing and memory usage for every individual job.


## 6. Advanced: SIDLE Multi-Region Reconstruction

### 6.1 When Does SIDLE Run?

SIDLE runs automatically when FAVABEAN detects **two or more amplicon regions** in your `files_info.csv`. You do not need to configure anything — the pipeline inspects the `region` column and, if it finds more than one unique value, triggers the SIDLE reconstruction after the per-region DADA2 processing completes.

If you have only one region, SIDLE is skipped entirely and the per-region ASV tables are your final output.

### 6.2 SIDLE Configuration

The `sidle:` section in `favabean.yaml` controls the reconstruction:

```yaml
sidle:
  max_mismatch: 2           # Mismatches allowed in ASV-to-reference alignment
  primer_max_mismatch: 1    # Mismatches allowed in in-silico PCR
  min_counts: 1000          # Minimum reads per sample for EM
  min_abund: 1.0e-5         # Minimum abundance threshold for convergence
  region_normalize:
    - average               # How to handle taxa found in multiple regions
```

**`region_normalize` options:**

- `average` (default) — Divides each feature's counts by the number of regions it was detected in. Prevents over-counting taxa present in multiple regions.
- `weighted` — Scales counts by (total regions / regions detected). Taxa found in fewer regions are scaled up.
- `unweighted` — No normalization. Raw EM-reconstructed counts are used as-is.

### 6.3 SIDLE Outputs

SIDLE produces **two pairs** of output tables, both in TSV and BIOM format:

**Accession-level** (`reconstructed_ASV.tsv`):

```
#OTU ID           P1A1    P1A2    taxonomy
018_0511          19.0    11.0    Bacteria;Proteobacteria;...;Neisseria sp._HMT_018
046_4327          133.0   128.0   Bacteria;Firmicutes;...;Gemella morbillorum
056AY020          125.0   152.0   Bacteria;Firmicutes;...;Streptococcus sp._HMT_056
058BM035          27.0    33.0    Bacteria;Firmicutes;...;Streptococcus oralis_subsp._dentisani_clade_058
```

Each row corresponds to a unique reference sequence (clade) in the database. The OTU IDs are accession numbers from the reference database. This preserves the full granularity of the SMURF EM reconstruction — if two reference sequences have different accession numbers but the same species name, they appear as separate rows.

**Species-condensed** (`reconstructed_ASV_species.tsv`):

```
#OTU ID                                                  P1A1    P1A2    taxonomy
Bacteria;Proteobacteria;...;Neisseria sp._HMT_018       19.0    11.0    ...
Bacteria;Firmicutes;...;Gemella morbillorum              133.0   128.0   ...
Bacteria;Firmicutes;...;Streptococcus sp._HMT_056       125.0   152.0   ...
```

This table merges rows that share the same full taxonomy string by summing their counts. The OTU ID is the taxonomy string itself. This is typically what you want for downstream analysis (alpha diversity, differential abundance testing, etc.) because it groups at the species level rather than the individual reference clade level.

In the example dataset, SIDLE starts with 80 unique reference clades (accession-level) and condenses them to 62 unique taxa (species-level).

### 6.4 SIDLE Summary Report

`preprocessing_summary_sidle.tsv` provides a concise overview of the reconstruction:

| taxonomy_database | n_regions_combined | regions | unique_taxa_after\_sidle\_reconstruct | pct_species_after\_sidle\_reconstruct | unique_taxa_after\_sidle\_to\_biom | pct_species_after\_sidle\_to\_biom |
|---|---|---|---|---|---|---|
| eHOMD_V15.22 | 2 | 27F,V3V5 | 80 | 100 | 62 | 100 |

- `unique_taxa_after_sidle_reconstruct` — Number of reference clades with non-zero counts after EM reconstruction (accession-level).
- `pct_species_after_sidle_reconstruct` — Percentage of those clades that have species-level taxonomy.
- `unique_taxa_after_sidle_to_biom` — Number of unique taxa after species-level condensation.
- `pct_species_after_sidle_to_biom` — Percentage of condensed taxa at species level.

The reduction from 80 to 62 means that 18 reference clades shared taxonomy strings with other clades and were merged during species condensation.

### 6.5 SIDLE Intermediate Files

The `favabean/sidle/` directory contains intermediate files from the reconstruction:

```
favabean/sidle/
├-- reference/
│   ├-- prepared_ref.fasta          # Reference FASTA with unique accession IDs
│   ├-- taxonomy.tsv                # Accession-to-taxonomy mapping
│   ├-- 27F_extracted.fasta         # In-silico PCR products for 27F region
│   ├-- 27F_kmers.fasta             # Fixed-length k-mers for alignment
│   ├-- 27F_kmer_map.tsv            # K-mer to parent sequence mapping
│   ├-- V3V5_extracted.fasta        # Same for V3V5 region
│   ├-- V3V5_kmers.fasta
│   └-- V3V5_kmer_map.tsv
├-- 27F/
│   ├-- rep_seqs.fasta              # ASV representative sequences
│   ├-- asv_counts.tsv              # ASV count matrix
│   ├-- seqlengths.txt              # Trim length for this region
│   ├-- trimmed_seqs.fasta          # Uniform-length ASVs
│   ├-- trimmed_counts.tsv          # Counts for trimmed ASVs
│   └-- alignment.tsv               # ASV-to-reference k-mer alignments
├-- V3V5/
│   └-- (same structure as 27F/)
├-- database_map.tsv                # Cross-region database map
├-- database_summary.tsv            # Database coverage summary
├-- reconstructed_counts.tsv        # Raw EM output (before taxonomy merge)
└-- reconstructed_taxonomy.tsv      # Taxonomy for reconstructed features
```


## 7. Using Outputs with FALAPhyl

FAVABEAN produces BIOM files that can be directly used with [FALAPhyl](https://github.com/khalidtab/FALAPhyl) for downstream statistical analysis (alpha diversity, beta diversity, differential abundance testing). Copy the BIOM files to your FALAPhyl `data/` directory:

```bash
cp favabean/*_ASV.biom /path/to/FALAPhyl/data/
# Or for the multi-region reconstructed table:
cp favabean/reconstructed_ASV_species.biom /path/to/FALAPhyl/data/
```

Then run FALAPhyl:

```bash
snakemake alpha beta diff --use-conda --cores all
```

The BIOM files are also compatible with other 16S analysis tools: QIIME2 (`qiime tools import`), phyloseq (R), microViz, MaAsLin2, etc.


## 8. Troubleshooting

**Very low chimera retention (<30%)**

Some chimera loss is normal, especially for longer amplicons (V1–V3 can lose 40–50% of reads to chimeras). If retention is extremely low, check that your `expected length` in `files_info.csv` is correct. An incorrect expected length causes Figaro to choose poor trim parameters, which leads to poor paired-end merging and high chimera rates.

**Conda environment creation fails**

If you are behind a firewall or proxy, Conda may fail to download packages. Try setting Conda to use your proxy, or use the Docker image which has all environments pre-built.

**SIDLE did not produce output**

SIDLE only runs when there are 2+ amplicon regions. Check that your `files_info.csv` has at least two distinct values in the `region` column. Also check `data/logs/sidle-*.log` for error messages.

**Pipeline interrupted mid-run**

Simply re-run the same `snakemake` command. Snakemake tracks all completed outputs and only re-runs the steps that did not finish. To force re-execution from a specific step, delete that step's output files. Sometimes after interrupting the pipeline, the folder will be "locked". To unlock the folder and rerun the pipeline, use `--unlock` in your command, then rerun it again without that flag.

**Out of memory**

Taxonomy assignment loads the full reference database into memory. eHOMD requires ~770 MB; SILVA requires ~4 GB. If your system is memory-constrained, reduce the `--cores` parameter to limit parallelism, which reduces peak memory usage.


## 9. FAQ

**Can I use databases other than eHOMD or SILVA?**

Yes. Any DADA2-formatted reference database works. You need two files: a taxonomy-training FASTA (taxonomy as header, used for genus-level classification) and a species-assignment FASTA (accession + binomial, used for species-level classification). Add a new entry under `taxonomy_database:` in `favabean.yaml` with the download URLs or local paths.

**What if I have only one amplicon region?**

The pipeline works identically, but SIDLE is not triggered. Your final outputs are the per-region ASV and taxonomy tables (e.g., `V4_eHOMD_ASV.tsv`, `V4_eHOMD_ASV.biom`). No `reconstructed_*` files are produced.

**Can I add more samples later?**

Yes. Add the new rows to `files_info.csv`, delete `data/files_info_Batches.csv` (so batch detection re-runs), and re-execute the pipeline. As long as you still have the files from the previous run, Snakemake will only process the new samples and regenerate the merged outputs. 

**What is the difference between `reconstructed_ASV` and `reconstructed_ASV_species`?**

`reconstructed_ASV` has one row per reference database clade (e.g., different strains or sub-species of the same organism appear as separate rows). `reconstructed_ASV_species` merges rows that share identical taxonomy strings by summing their counts. Use the species-condensed table for most downstream analyses; use the accession-level table if you need strain-level resolution, or if your analysis requires you using the phylogenic tree produced by the taxonomy database.

**Why do reads_after_sidle_reconstruct differ from per-region counts?**

SIDLE's EM algorithm does not simply sum the per-region reads. It probabilistically assigns ambiguous reads to reference organisms using evidence from all regions simultaneously. This redistribution means the total can be higher or lower than any single region's total. This is expected behavior and reflects the improved taxonomic resolution from combining multiple regions. That being said, if the same sequence was assigned the same taxonomic rank across multiple regions, the default behaviour of the SIDLE algorithm, as referenced in `favabean'yaml` is to average the two counts `region_normalize: average`.

**How do I interpret the `Batch_ID` in files_info_Batches.csv?**

This file is auto-generated by `batch.R`, which reads the FASTQ headers to detect sequencing batches. Each unique combination of instrument ID and flowcell ID becomes a batch. Instructions on how to invoke this script can be found when you write the command `start`, which will show you how to initialize and create your `files_info_Batches.csv`  file.  Meaning, you never need to edit this file manually, and need to know which batch each sequence file comes from.
