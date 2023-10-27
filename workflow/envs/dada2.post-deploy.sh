#!/bin/bash

# Initialize conda
. /opt/conda/etc/profile.d/conda.sh

# Activate the specific conda environment
conda activate $CONDA_PREFIX

# Navigate to the desired directory
cd ./workflow/envs/

# Run the R script to install packages
Rscript install_packages.R
