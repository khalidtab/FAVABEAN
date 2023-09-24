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
        "data/cutadapt/{batch}-{region}/.{sample}-cutadapt.done"
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
        mkdir -p /data/cutadapt/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/ &&
        mkdir -p /data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/ &&
        cutadapt -j 0 --minimum-length {params.filt} -g {params.primer_5} -G {params.primer_3} \
        --too-short-output        /data/cutadapt/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq \
        --too-short-paired-output /data/cutadapt/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq \
        -o /data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq \
        -p /data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq \
        {input.fq1} {input.fq2} > {log} 2>&1 && touch {output}
        """

rule seqkit:
    input:
        lambda wildcards: expand("data/cutadapt/{batch}-{region}/.{sample}-cutadapt.done",
                                 batch=wildcards.batch,
                                 region=wildcards.region,
                                 sample=get_samples_from_batch_region(wildcards.batch, wildcards.region))
    output:
        R1="data/cutadapt/{batch}-{region}/stats_R1_lengths.done",
        R2="data/cutadapt/{batch}-{region}/stats_R2_lengths.done",
        R1trim=touch("data/cutadapt/{batch}-{region}/.trim_R1.done"),
        R2trim=touch("data/cutadapt/{batch}-{region}/.trim_R2.done")
    log:
        R1="data/logs/seqkit_BBMAP-{batch}-{region}R1.log",
        R2="data/logs/seqkit_BBMAP-{batch}-{region}R2.log"
    params:
        filt=config["initial_filter"],
        trim=config["trim_param"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # R1
        # If its a number, this will be placed here. If it is not, it will be rewritten in one of the codes below
        echo {params.trim} > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1{params.trim}.txt
        cat data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/*R1_001.fastq | seqkit stats > {output.R1} -T -a && 
        values=$(tail -n 1 {output.R1} | cut -f7-11)
        echo $values | cut -d' ' -f2 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1max_len.txt
        echo $values | cut -d' ' -f3 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1Q1.txt
        echo $values | cut -d' ' -f4 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1Q2.txt
        echo $values | cut -d' ' -f5 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1Q3.txt
        
        x_int=$(printf "%.0f\n" $(cat data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1Q1.txt))
        y_int=$(printf "%.0f\n" $(cat data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1Q2.txt))
        ten_percent=$((y_int / 10))
        lower_bound=$((y_int - ten_percent))
        upper_bound=$((y_int + ten_percent))
        if [ "$x_int" -ge "$lower_bound" ] && [ "$x_int" -le "$upper_bound" ]; then
            echo $x_int > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1default.txt
        else
            echo $y_int > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1default.txt
        fi
        
        trim_file=$(cat "data/cutadapt/{wildcards.batch}-{wildcards.region}/.R1{params.trim}.txt")
        
        ls /data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/*R1*.fastq | parallel "mkdir -p {{//}}/trim/ && bbduk.sh in={{}} minlength=$trim_file maxlength=$trim_file nullifybrokenquality out=stdout.fq 2>> {log.R1} | seqkit seq --remove-gaps -m $trim_file -M $trim_file > {{//}}/trim/{{/}}"
        
        # R2
        # If its a number, this will be placed here. If it is not, it will be rewritten in one of the codes below
        echo {params.trim} > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2{params.trim}.txt
        cat data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/*R2_001.fastq | seqkit stats > {output.R2} -T -a && 
        values=$(tail -n 1 {output.R2} | cut -f7-11)
        echo $values | cut -d' ' -f2 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2max_len.txt
        echo $values | cut -d' ' -f3 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2Q1.txt
        echo $values | cut -d' ' -f4 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2Q2.txt
        echo $values | cut -d' ' -f5 > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2Q3.txt
        x_int=$(printf "%.0f\n" $(cat data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2Q1.txt))
        y_int=$(printf "%.0f\n" $(cat data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2Q2.txt))
        ten_percent=$((y_int / 10))
        lower_bound=$((y_int - ten_percent))
        upper_bound=$((y_int + ten_percent))
        if [ "$x_int" -ge "$lower_bound" ] && [ "$x_int" -le "$upper_bound" ]; then
            echo $x_int > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2default.txt
        else
            echo $y_int > data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2default.txt
        fi
        
        trim_file=$(cat "data/cutadapt/{wildcards.batch}-{wildcards.region}/.R2{params.trim}.txt")
        
        ls data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/*R2*.fastq | parallel "mkdir -p {{//}}/trim/ && bbduk.sh in={{}} minlength=$trim_file maxlength=$trim_file nullifybrokenquality out=stdout.fq 2>> {log.R2} | seqkit seq --remove-gaps -m $trim_file -M $trim_file > {{//}}/trim/{{/}}"
        
        """

rule figaro:
    input:
        "data/cutadapt/{batch}-{region}/.trim_R2.done"
    output:
        figaroFold=directory("data/cutadapt/{batch}-{region}/figaro"),
        figaroDone="data/cutadapt/{batch}-{region}/.figaro.txt"
    log:
        R1="data/logs/figaro-{batch}-{region}.log"
    conda:
        "../envs/figaro.yaml"
    shell:
        """
         python workflow/envs/figaro/figaro/figaro.py -i data/cutadapt/{wildcards.batch}-{wildcards.region}/cutadapt/trim -o {output.figaroFold} -f 1 -r 1 -a 400 -F illumina > {log}
         touch {output.figaroDone}
        """




combinations = [(row['Batch_ID'], row['region']) for _, row in samples_table.drop_duplicates(['Batch_ID', 'region']).iterrows()]

rule all:
    input:
        expand("data/cutadapt/{batch}-{region}/.figaro.txt", 
                batch=[combo[0] for combo in combinations],
               region=[combo[1] for combo in combinations])

