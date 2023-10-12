import pandas as pd
# Read the CSV file
samples_table = pd.read_csv("data/files_info.csv").set_index(["sample", "region"], drop=False)

# Internally adjust the `fastq1` and `fastq2` columns
samples_table["fastq1"] = "data/" + samples_table["fastq1"]
samples_table["fastq2"] = "data/" + samples_table["fastq2"]

def check_column_names(df):
    columns = df.columns
    stripped_columns = [col.strip() for col in columns]
    whitespace_columns = [col for col, stripped_col in zip(columns, stripped_columns) if col != stripped_col]
    return whitespace_columns

whitespace_columns = check_column_names(samples_table)
if whitespace_columns:
    print(f"Warning: The following columns have leading or trailing whitespaces: {whitespace_columns}")

def get_samples_from_batch_region(batch, region):
    return samples_table.loc[(samples_table['Batch_ID'] == batch) & (samples_table['region'] == region), 'sample'].tolist()

class Extractor:
    def __init__(self):
        pass

    def getSampleInfo(self, fileName:str):
        baseName = fileName.split(".")[0]
        baseSplit = baseName.split("_")
        sampleNum = int(baseSplit[-4].replace("S",""))
        return sampleNum

# Create an instance of the Extractor class
extractor = Extractor()

# Apply the getSampleInfo method to the fastq1 column and store the results in a new column
samples_table["SampleNum"] = samples_table["fastq1"].apply(lambda x: extractor.getSampleInfo(x))

# Now the samples_table variable has the updated "SampleNum" column and can be used in the rest of your script.

rule cutadapt:
    input:
        fq1 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq1'],
        fq2 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq2']
    output:
        touch("data/favabean/{batch}-{region}/.{sample}-cutadapt.done")
    log:
        "data/logs/cutadapt-{batch}-{region}-{sample}.log"
    params:
        primer_5 =      lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_5'],
        primer_3 =      lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_3'],
        sampleNum =     lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'SampleNum'],
        filt=config["initial_filter"]
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        mkdir -p /data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/ &&
        mkdir -p /data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/ &&
        cutadapt -j 0 --minimum-length {params.filt} -g {params.primer_5} -G {params.primer_3} \
        --too-short-output        /data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq \
        --too-short-paired-output /data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq \
        -o /data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq \
        -p /data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq \
        {input.fq1} {input.fq2} > {log} 2>&1 
        """





rule seqkitR1:
    input:
        lambda wildcards: expand("data/favabean/{batch}-{region}/.{sample}-cutadapt.done",
                                 batch=wildcards.batch,
                                 region=wildcards.region,
                                 sample=get_samples_from_batch_region(wildcards.batch, wildcards.region))
    output:
        R1done=touch("data/favabean/{batch}-{region}/.trim_R1.done")
    log:
        R1="data/logs/seqkit_BBMAP-{batch}-{region}R1.log"
    params:
        filt=config["initial_filter"],
        trim=config["trim_param"],
        R1="R1"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # R1
        # If its a number, this will be placed here. If it is not, it will be rewritten in one of the codes below
        echo {params.trim} > /data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/.{params.R1}{params.trim}.txt
        cat data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/*{params.R1}_001.fastq | seqkit stats >> data/favabean/{wildcards.batch}-{wildcards.region}/stats_{params.R1}_lengths.txt -T -a && 
        values=$(tail -n 1 data/favabean/{wildcards.batch}-{wildcards.region}/stats_{params.R1}_lengths.txt | cut -f7-11)
        echo $values | cut -d' ' -f2 > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}max_len.txt
        echo $values | cut -d' ' -f3 > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q1.txt
        echo $values | cut -d' ' -f4 > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q2.txt
        echo $values | cut -d' ' -f5 > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q3.txt
        
        x_int=$(printf "%.0f\n" $(cat /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q1.txt))
        y_int=$(printf "%.0f\n" $(cat /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q2.txt))
        ten_percent=$((y_int / 10))
        lower_bound=$((y_int - ten_percent))
        upper_bound=$((y_int + ten_percent))
        if [ "$x_int" -ge "$lower_bound" ] && [ "$x_int" -le "$upper_bound" ]; then
            echo $x_int > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}default.txt
        else
            echo $y_int > /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}default.txt
        fi
        
        trim_file=$(cat "/data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}{params.trim}.txt")
        
        ls /data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/*{params.R1}*.fastq | parallel "mkdir -p {{//}}/trim/ && bbduk.sh in={{}} minlength=$trim_file maxlength=$trim_file nullifybrokenquality out=stdout.fq 2>> {log.R1} | seqkit seq --remove-gaps -m $trim_file -M $trim_file > {{//}}/trim/{{/}}"
        rm -f /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}{params.trim}.txt /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}max_len.txt /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q1.txt /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q2.txt /data/favabean/{wildcards.batch}-{wildcards.region}/.{params.R1}Q3.txt
        """



use rule seqkitR1 as seqkitR2 with:
    output:
        R1done=touch("data/favabean/{batch}-{region}/.trim_R2.done")
    log:
        R1="data/logs/seqkit_BBMAP-{batch}-{region}R2.log"
    params:
        filt=config["initial_filter"],
        trim=config["trim_param"],
        R1="R2"



#    return min(2 * core_multiplier, total_cores)



rule figaro:
    input:
        "data/favabean/{batch}-{region}/.trim_R1.done",
        "data/favabean/{batch}-{region}/.trim_R2.done"
    output:
        "data/favabean/{batch}-{region}/figaro/trimParameters.json"
    log:
        R1="data/logs/figaro-{batch}-{region}.log"
    conda:
        "../envs/figaro.yaml"
    shell:
        """
         python workflow/envs/figaro/figaro/figaro.py -i data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/trim -o data/favabean/{wildcards.batch}-{wildcards.region}/figaro/ -f 1 -r 1 -a 500 -F illumina >> {log}
         # Do some clean up of files
         rm -rf data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/trim
        """



def determine_threads(wildcards):
    total_cores = workflow.cores
    core_multiplier = workflow.attempt
    return total_cores
    # This sets a base of 2 threads per job and adjusts based on the attempt number.
    # Adjust this logic as needed for your specific use-case.


# Note that this uses the cutadapt files, and not the figaro ones!
rule dada2:
    input:
        "data/favabean/{batch}-{region}/figaro/trimParameters.json"
    output:
        "data/favabean/{batch}-{region}/cutadapt/seqtab.tsv"
    log:
        R1="data/logs/dada2-{batch}-{region}.log"
    conda:
        "../envs/dada2.yaml"
    threads:
        determine_threads
    params:
        figaro=config["figaro"]
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2.R -i {input} -f data/favabean/{wildcards.batch}-{wildcards.region}/ -p {params.figaro} -c {threads} >> {log} 2>&1
        """


combinations = [(row['Batch_ID'], row['region']) for _, row in samples_table.drop_duplicates(['Batch_ID', 'region']).iterrows()]

rule all:
    input:
        expand("data/favabean/{batch}-{region}/cutadapt/seqtab.tsv", 
                batch=[combo[0] for combo in combinations],
               region=[combo[1] for combo in combinations])

