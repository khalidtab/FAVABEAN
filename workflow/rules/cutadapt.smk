import pandas as pd

# Read the CSV file
samples_table = pd.read_csv("data/files_info.csv").set_index(["sample", "region"], drop=False)

# Internally adjust the `fastq1` and `fastq2` columns
samples_table["fastq1"] = "data/" + samples_table["fastq1"]
samples_table["fastq2"] = "data/" + samples_table["fastq2"]

# Check for leading or trailing whitespaces in column names
import re

def check_column_names(df):
    columns = df.columns
    stripped_columns = [col.strip() for col in columns]
    whitespace_columns = [col for col, stripped_col in zip(columns, stripped_columns) if col != stripped_col]

    return whitespace_columns

whitespace_columns = check_column_names(samples_table)

if whitespace_columns:
    print(f"Warning: The following columns have leading or trailing whitespaces: {whitespace_columns}")

# Function to generate input files for cutadapt rule
def generate_cutadapt_inputs():
    inputs = []
    for _, row in samples_table.iterrows():
        inputs.extend([
            f"data/cutadapt/{row['Batch_ID']}-{row['region']}/{row['sample']}--R1.fastq",
            f"data/cutadapt/{row['Batch_ID']}-{row['region']}/{row['sample']}--R2.fastq"
        ])
    return inputs

rule cutadapt:
    input:
        lambda wildcards: generate_cutadapt_inputs()

rule cutadapt_action:
    input:
        fq1 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq1'],
        fq2 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq2']
    output:
        R1 = "data/cutadapt/{batch}-{region}/{sample}--R1.fastq",
        R2 = "data/cutadapt/{batch}-{region}/{sample}--R2.fastq"
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
