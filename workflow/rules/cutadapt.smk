import pandas as pd

samples_table = pd.read_csv("data/files_info.csv").set_index(["sample", "region"], drop=False)

# Add this code after reading the CSV file
import re

# Check for leading or trailing whitespaces in column names
def check_column_names(df):
    columns = df.columns
    stripped_columns = [col.strip() for col in columns]
    whitespace_columns = [col for col, stripped_col in zip(columns, stripped_columns) if col != stripped_col]

    return whitespace_columns

whitespace_columns = check_column_names(samples_table)

if whitespace_columns:
    print(f"Warning: The following columns have leading or trailing whitespaces: {whitespace_columns}")

# Rest of your Snakemake code...




# Function to return fq1 and fq2 for a given sample-region combination.
def fq_dict_from_sample(wildcards):
    return {
        "fq1": samples_table.loc[(wildcards.sample, wildcards.region), "fastq1"],
        "fq2": samples_table.loc[(wildcards.sample, wildcards.region), "fastq2"]
    }

rule cutadapt:
    input:
        expand(["data/cutadapt/cutadapt-{batch}-{sample}-{region}-R1.fastq",
                "data/cutadapt/cutadapt-{batch}-{sample}-{region}-R2.fastq"],
               batch=samples_table["Batch_ID"].unique(), 
               sample=samples_table["sample"].unique(), 
               region=samples_table["region"].unique())

rule cutadapt_action:
    input:
        fq1 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq1'],
        fq2 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq2']
    output:
        R1 = "data/cutadapt/cutadapt-{batch}-{sample}-{region}-R1.fastq",
        R2 = "data/cutadapt/cutadapt-{batch}-{sample}-{region}-R2.fastq"
    log:
        "data/logs/cutadapt-{batch}-{sample}-{region}.log"
    params:
        primer_5 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_5'],
        primer_3 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_3']
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        cutadapt -g {params.primer_5} -G {params.primer_3} -o {output.R1} -p {output.R2} {input.fq1} {input.fq2} --minimum-length 50 > {log} 2>&1
        """
