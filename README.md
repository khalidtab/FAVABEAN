[![GitHub Super-Linter](https://github.com/khalidtab/FALAPhyl/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)

<div style="text-align: center;">
  <img src="FAVABEAN.png" alt="Banner" style="max-width: 100%; height: auto;">
</div>

# FAVABEAN: *F*ast *A*mplicon *V*ariant *A*nnotation, *B*inning, *E*rror-correction and *AN*alysis
This is a pipeline that fully automates creating biom files from fastq files, by running Figaro to find the best parameters for merging in DADA2, applying DADA2, then running a naive Bayesian taxonomy identification.

## Example
Pipeline is best used with its Docker container, using the following command:

> docker run --rm -it -v ~/Path/To/Your/Folder:/data khalidtab/favabean:latest start

This will copy the following files to your path (if they are not there already):

- `FAVABEAN_environments.txt`: Contains the paths to the snakemake environments used by the pipeline. It is helpful if you want to run specific commands on the environment yourself.
- `files_info.csv`: template file to use
- `favabean.yaml`: configuration file

If this is the first time you are running the pipeline on this dataset, you will need to initialize your `files_info.csv` file. This is done by the following:

> conda activate <<biom environment such as: .snakemake/conda/007a7beaa3353a33c938a5a0e57be4ff_, this can be found in the FAVABEAN_environments.txt file>>

Next, run the following command:
> Rscript --vanilla workflow/scripts/batch.R

The `batch.R` file will read all your fastq/gz files you have listed in your `files_info.csv`, and designate them into their own separate sequencing batches based on the information written in the files. If successful, you will have a file called `files_info_Batches.csv` generated in your folder, with all that information.

If all is ready, you can run the pipeline by executing the following command (which, is conveniently printed into your console as well)

> snakemake paired_taxonomy --use-conda --cores all --keep-going --retries 5 --rerun-incomplete --scheduler greedy

Note that `--use-conda` is mandatory for the correct execution of the pipeline, but the rest of the commands are modifiable. Check `snakemake --help` for more information.

