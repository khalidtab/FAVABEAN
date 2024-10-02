#!/bin/bash

# Step 1: Copy files as needed
cp -n ./files_info.csv data/
cp -n ./input.yaml data/ 
cp -n ./environments.txt data/ 


# Step 2: Pre-populate the prompt with the Snakemake command
# Use 'bind' to set the command in the shell's buffer
bind 'set pre-input-hook "echo snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete; bind -r pre-input-hook"'

# Step 3: Start bash
exec bash