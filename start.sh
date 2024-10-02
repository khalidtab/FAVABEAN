#!/bin/bash

echo "Welcome to FAVABEAN!"
echo ""
echo ""
echo "Copying the template files in your directly if they don't exist (environments.txt, files_info.csv, and input.yaml)'"
# Step 1: Copy files as needed (if applicable)
cp -p environments.txt data/
cp -p files_info.csv data/
cp -p input.yaml data/

echo ""
echo ""
echo "The following command(s) can be executed:"
echo "snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete"
echo "Please copy and paste the above command or type your own command below."

# Open an interactive bash shell so the user can input commands
exec bash