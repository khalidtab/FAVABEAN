#!/bin/bash
echo ""
echo ""
echo -e "\033[0;32m====================\033[0m"
echo -e "\033[0;32mWelcome to FAVABEAN!\033[0m"
echo -e "\033[0;32m====================\033[0m"
echo ""
echo ""
echo "Copying the template files in your directly if they don't exist (environments.txt, files_info.csv, and favabean.yaml)'"
# Step 1: Copy files as needed (if applicable)
cp -pn environments.txt data/FAVABEAN_environments.txt
cp -pn files_info.csv data/
cp -pn favabean.yaml data/

echo "If this is the first time you are running the pipeline on this dataset, you will need to initialize your 'files_info.csv' file."
echo "Please write 'conda activate' and copy the biom.yaml environment link from FAVABEAN_environments.txt file"
echo ""
echo "Example:"
echo -e "\033[0;36m conda activate <<biom environment such as : .snakemake/conda/007a7beaa3353a33c938a5a0e57be4ff_>> \033[0m"
echo ""
echo "Then, execute the following command"
echo -e "\033[0;36msRscript --vanilla workflow/scripts/batch.R\033[0m"
echo ""
echo ""
echo "Example command(s):"
echo -e "\033[0;36msnakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete --scheduler greedy \033[0m"
echo "Please copy and paste the above command(s) or type your own command below."
echo ""
echo ""

# Open an interactive bash shell so the user can input commands
exec bash
