#!/bin/bash

# Step 1: Copy files as needed (if applicable)
cp /path/to/source /path/to/destination

# Step 2: Display the Snakemake command but don't execute it
echo "The following command is ready for execution:"
echo "snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete"

# Step 3: Prompt the user to press enter to pre-fill the command
echo "Press Enter to execute the Snakemake command or type your own command below:"
read -e -p "snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete" user_input

# If the user presses enter without typing anything, execute the default command
if [ -z "$user_input" ]; then
    snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete
else
    # Execute the user-specified command
    $user_input
fi

# Step 4: Start bash for further interaction (optional, if you want to keep the shell open)
exec bash