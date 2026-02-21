#!/bin/bash
set -euo pipefail

# =============================================================================
# FAVABEAN: Automated Working Example
# =============================================================================
# This script sets up and runs the FAVABEAN pipeline end-to-end using
# small bundled FASTQ test files. No internet connection is required
# for the FASTQ data — only for Conda environment setup and taxonomy
# database downloads (handled automatically by the pipeline).
#
# Test dataset:
#   - 2 biological samples (P1A1, P1A2)
#   - 2 amplicon regions (27F targeting V1-V3, V3V5 targeting V3-V5)
#   - 2 sequencing runs (different Illumina flowcells)
#   - 16 FASTQ files total (2000 reads each, ~1.7 MB)
#   This demonstrates automatic batch detection across instrument runs,
#   multi-region processing, and per-region taxonomy assignment.
#
# Usage:
#   cd /path/to/FAVABEAN
#   bash examples/run_example.sh
#
# Prerequisites:
#   - Snakemake (>= 7.0)
#   - Conda or Mamba
#
# What this script does:
#   1. Copies example FASTQ files and config into data/
#   2. Runs batch detection (determines sequencing batches from headers)
#   3. Executes the full FAVABEAN pipeline via Snakemake
#   4. Reports the output location
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "============================================="
echo " FAVABEAN: Working Example"
echo "============================================="
echo ""
echo " Project directory: $PROJECT_DIR"
echo ""

# ---- Step 0: Verify we're in the right place ----
if [[ ! -f "$PROJECT_DIR/workflow/snakefile" ]]; then
    echo "ERROR: Cannot find workflow/snakefile."
    echo "Please run this script from the FAVABEAN root directory:"
    echo "  bash examples/run_example.sh"
    exit 1
fi

# ---- Step 1: Set up data directory ----
echo "Step 1: Setting up data directory..."

# Safety check — don't overwrite existing data
if [[ -d "$PROJECT_DIR/data" ]]; then
    echo ""
    echo "WARNING: data/ directory already exists."
    echo "This example will overwrite files_info.csv and favabean.yaml."
    read -p "Continue? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
fi

mkdir -p "$PROJECT_DIR/data"

# Copy example FASTQ files
echo "  Copying example FASTQ files to data/..."
cp "$SCRIPT_DIR"/test_data/*.fastq.gz "$PROJECT_DIR/data/"

# Copy metadata
echo "  Copying files_info.csv..."
cp "$SCRIPT_DIR/files_info.csv" "$PROJECT_DIR/data/files_info.csv"

# Copy config
echo "  Copying favabean.yaml..."
cp "$SCRIPT_DIR/favabean.yaml" "$PROJECT_DIR/data/favabean.yaml"

echo "  Done."
echo ""

# ---- Step 2: Run the pipeline ----
echo "Step 2: Running FAVABEAN pipeline..."
echo "  (This will create Conda environments on first run — may take a while)"
echo ""

cd "$PROJECT_DIR"

# Remove any stale batch file so it gets regenerated
rm -f data/files_info_Batches.csv

snakemake paired \
    --snakefile workflow/snakefile \
    --use-conda \
    --cores all \
    --keep-going \
    --retries 2 \
    --rerun-incomplete \
    --scheduler greedy

echo ""
echo "============================================="
echo " FAVABEAN example complete!"
echo "============================================="
echo ""
echo " Output files are in: data/favabean/"
echo ""
echo " Key outputs:"
echo "   - *_OTU.tsv           — OTU/ASV count tables"
echo "   - *_taxonomy.tsv      — Taxonomy assignments"
echo "   - *.biom              — BIOM-format tables"
echo "   - preprocessing_summary.tsv — Pipeline summary statistics"
echo ""
echo " To run downstream analysis with FALAPhyl, copy the"
echo " BIOM files to your FALAPhyl data/ directory and run:"
echo "   snakemake alpha beta diff --use-conda --cores all"
echo "============================================="
