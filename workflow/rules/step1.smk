import subprocess
import os
import pandas as pd

# Path to the output CSV file generated by the R script
output_csv = "data/files_info_Batches.csv"

# Step 1: Check if the file already exists
if not os.path.exists(output_csv):
    try:
        # Run the R script if the file does not exist
        print("Sequencing batches information not available. Reading the qzipped files and determining that now.")
        subprocess.run(
            f"conda run -p /.snakemake/conda/ee5771f5a5b4f769c2b98bafaa645bd0_ Rscript workflow/scripts/batch.R",
            shell=True, check=True)
        print("Sequencing batches determined. File saved as 'files_info_Batches.csv'")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running the R script: {e}")
        exit(1)
else:
    print(f"Sequencing batch  information already exists, as the {output_csv} already exists,  skipping R script execution.")

# Step 2: Read the CSV file (whether generated or already existing)
samples_table = pd.read_csv(output_csv).set_index(["sample", "region"], drop=False)

print(f"Loaded samples table with {len(samples_table)} entries.")


# Internally adjust the `fastq1` and `fastq2` columns
samples_table["fastq1"] = "data/" + samples_table["fastq1"]
samples_table["fastq2"] = "data/" + samples_table["fastq2"]

# Check for column names with leading or trailing whitespaces
def check_column_names(df):
    columns = df.columns
    stripped_columns = [col.strip() for col in columns]
    whitespace_columns = [col for col, stripped_col in zip(columns, stripped_columns) if col != stripped_col]
    return whitespace_columns

whitespace_columns = check_column_names(samples_table)
if whitespace_columns:
    print(f"Warning: The following columns have leading or trailing whitespaces: {whitespace_columns}")

# Get samples from batch and region
def get_samples_from_batch_region(batch, region):
    return samples_table.loc[(samples_table['Batch_ID'] == batch) & (samples_table['region'] == region), 'sample'].tolist()

# Get expected length from batch and region
def get_expected_length_from_batch_region(batch, region):
    lengths = samples_table.loc[
        (samples_table['Batch_ID'] == batch) & (samples_table['region'] == region), 'expected.length'
    ].unique()
    if len(lengths) == 1:
        return lengths[0]
    else:
        raise ValueError(f"Multiple expected lengths found for batch {batch} and region {region}.")


# Define thread management functions
def determine_threads(wildcards):
    total_cores = workflow.cores
    core_multiplier = workflow.attempt
    return total_cores / 2

def all_threads(wildcards):
    total_cores = workflow.cores
    core_multiplier = workflow.attempt
    return total_cores
    
    
    #    return min(2 * core_multiplier, total_cores)
    # This sets a base of 2 threads per job and adjusts based on the attempt number.
    # Adjust this logic as needed for your specific use-case.
    
# Now the samples_table variable has the updated "SampleNum" column and can be used in the rest of your script.

rule cutadapt:
    input:
        fq1 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq1'],
        fq2 = lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'fastq2']
    output:
        touch("data/favabean/{batch}-{region}/.tmp/.{sample}-cutadapt.done")
    log:
        "data/logs/cutadapt-{batch}-{region}-{sample}.log"
    params:
        primer_5 =      lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_5'],
        primer_3 =      lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'primer_3'],
        sampleNum =     lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'SampleNum'],
        filt      =     config["initial_filter"]
    message: "Cutadapt - Removing the adaptors for sample {wildcards.sample} from batch: {wildcards.batch}, region: {wildcards.region}"
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        primer5=$(echo '{params.primer_5}' | sed 's/;/ -g /g')
        primer3=$(echo '{params.primer_3}' | sed 's/;/ -G /g')
        mkdir -p data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/ &&
        mkdir -p data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/ &&
        cutadapt -j 0 --minimum-length {params.filt} -g $primer5 -G $primer3 \
        --too-short-output        data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq.gz \
        --too-short-paired-output data/favabean/{wildcards.batch}-{wildcards.region}/initialFilt_discarded/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq.gz \
        -o data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq.gz \
        -p data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq.gz \
        {input.fq1} {input.fq2} > {log} 2>&1 
        """

rule sequence_length_stats:
    input:
        lambda wildcards: expand("data/favabean/{batch}-{region}/.tmp/.{sample}-cutadapt.done",
                                 batch=wildcards.batch,
                                 region=wildcards.region,
                                 sample=get_samples_from_batch_region(wildcards.batch, wildcards.region))
    output:
        R1="data/favabean/{batch}-{region}/.trimParamR1.txt",
        R2="data/favabean/{batch}-{region}/.trimParamR2.txt"
    params:
        filt=config["initial_filter"],
        trim_param=config["trim_param"]
    message: "SeqKit - Calculate descriptive statistics of the R1 and R2 sequence lengths for the samples: {wildcards.batch}, region: {wildcards.region}"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
            #R1
            echo {params.trim_param} > data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/.R1{params.trim_param}.txt
            cat data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/*R1_001.fastq.gz | gzip -d | seqkit stats > data/favabean/{wildcards.batch}-{wildcards.region}/stats_R1_lengths.txt -T -a
            
            values=$(tail -n 1 data/favabean/{wildcards.batch}-{wildcards.region}/stats_R1_lengths.txt | cut -f7-11)
            echo $values | cut -d' ' -f2 > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_max_len.txt
            echo $values | cut -d' ' -f3 > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_Q1.txt
            echo $values | cut -d' ' -f4 > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_Q2.txt
            echo $values | cut -d' ' -f5 > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_Q3.txt
            
            # Additional processing to determine default value
            x_int=$(printf "%.0f\n" $(cat data/favabean/{wildcards.batch}-{wildcards.region}/.R1_Q1.txt))
            y_int=$(printf "%.0f\n" $(cat data/favabean/{wildcards.batch}-{wildcards.region}/.R1_Q2.txt))
            ten_percent=$((y_int / 10))
            lower_bound=$((y_int - ten_percent))
            upper_bound=$((y_int + ten_percent))
            if [ "$x_int" -ge "$lower_bound" ] && [ "$x_int" -le "$upper_bound" ]; then
                echo $x_int > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_default.txt
            else
                echo $y_int > data/favabean/{wildcards.batch}-{wildcards.region}/.R1_default.txt
            fi
            cp data/favabean/{wildcards.batch}-{wildcards.region}/.R1_{params.trim_param}.txt data/favabean/{wildcards.batch}-{wildcards.region}/.trimParamR1.txt
            
            #R2
            echo {params.trim_param} > data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/.R2_{params.trim_param}.txt
            cat data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/*R2_001.fastq.gz | gzip -d | seqkit stats > data/favabean/{wildcards.batch}-{wildcards.region}/stats_R2_lengths.txt -T -a
            
            values=$(tail -n 1 data/favabean/{wildcards.batch}-{wildcards.region}/stats_R2_lengths.txt | cut -f7-11)
            echo $values | cut -d' ' -f2 > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_max_len.txt
            echo $values | cut -d' ' -f3 > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_Q1.txt
            echo $values | cut -d' ' -f4 > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_Q2.txt
            echo $values | cut -d' ' -f5 > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_Q3.txt
            
            # Additional processing to determine default value
            x_int=$(printf "%.0f\n" $(cat data/favabean/{wildcards.batch}-{wildcards.region}/.R2_Q1.txt))
            y_int=$(printf "%.0f\n" $(cat data/favabean/{wildcards.batch}-{wildcards.region}/.R2_Q2.txt))
            ten_percent=$((y_int / 10))
            lower_bound=$((y_int - ten_percent))
            upper_bound=$((y_int + ten_percent))
            if [ "$x_int" -ge "$lower_bound" ] && [ "$x_int" -le "$upper_bound" ]; then
                echo $x_int > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_default.txt
            else
                echo $y_int > data/favabean/{wildcards.batch}-{wildcards.region}/.R2_default.txt
            fi
            cp data/favabean/{wildcards.batch}-{wildcards.region}/.R2_{params.trim_param}.txt data/favabean/{wildcards.batch}-{wildcards.region}/.trimParamR2.txt

        """

rule cutAndKeepSameLengthSequencesForFigaro:
    input:
        sampleCutadapt=lambda wildcards: expand("data/favabean/{batch}-{region}/.tmp/.{sample}-cutadapt.done",
                                 batch=wildcards.batch,
                                 region=wildcards.region,
                                 sample=get_samples_from_batch_region(wildcards.batch, wildcards.region)),
        trimParamR1=lambda wildcards: expand("data/favabean/{batch}-{region}/.trimParamR1.txt",
                                 batch=wildcards.batch,
                                 region=wildcards.region),
        trimParamR2=lambda wildcards: expand("data/favabean/{batch}-{region}/.trimParamR2.txt",
                                 batch=wildcards.batch,
                                 region=wildcards.region)
    output:
        touch("data/favabean/{batch}-{region}/.{sample}_seqkit.done")
    log:
        "data/logs/seqkit-{batch}-{region}-{sample}.log"
    message: "SeqKit - Filtering Sample: {wildcards.sample} FASTQ files based on the chosen length for R1 and R2 from batch: {wildcards.batch}, region: {wildcards.region}"
    conda:
        "../envs/seqkit.yaml"
    params:
        sampleNum =     lambda wildcards: samples_table.loc[(wildcards.sample, wildcards.region), 'SampleNum']
    shell:
        """
        mkdir -p data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro/{wildcards.sample}
        
        trimParamR1=$(cat {input.trimParamR1})
        trimParamR2=$(cat {input.trimParamR2})
        
        # Cut sequences longer than the chosen trim parameter length, then remove gaps and any sequences shorter than that, then pair the sequences so that only those that are present in both would be kept
        cat data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_001.fastq.gz | gzip -d | seqkit subseq -r 1:$trimParamR1 | seqkit seq --remove-gaps -m $trimParamR1 -M $trimParamR1 > data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_trimmed.fastq
        cat data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_001.fastq.gz | gzip -d | seqkit subseq -r 1:$trimParamR2 | seqkit seq --remove-gaps -m $trimParamR2 -M $trimParamR2 > data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_trimmed.fastq
        seqkit pair -1 data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_trimmed.fastq -2 data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_trimmed.fastq -O data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro/{wildcards.sample} --force > {log} 2>&1
        # I think what seqkit pair kept on rewriting the whole folder at every execution, hence why I have each sample write their result in a separate folder. Now let's move its contents output
        mv data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro/{wildcards.sample}/* data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro/
        rm -r data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro/{wildcards.sample}/
        rm data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R1_trimmed.fastq data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/{wildcards.sample}_S{params.sampleNum}_L001_R2_trimmed.fastq
        """


rule figaro:
    input:
        lambda wildcards: [
            f"data/favabean/{wildcards.batch}-{wildcards.region}/.{sample}_seqkit.done"
            for sample in get_samples_from_batch_region(wildcards.batch, wildcards.region)
        ]
    output:
        "data/favabean/{batch}-{region}/figaro/trimParameters.json"
    log:
        R1="data/logs/figaro-{batch}-{region}.log"
    conda:
        "../envs/figaro.yaml"
    message: "Figaro - calculating the optimal parameters for DADA2 trimming for batch: {wildcards.batch}, region: {wildcards.region}"
    params:
        expected_length = lambda wildcards: get_expected_length_from_batch_region(wildcards.batch, wildcards.region)
    threads:
        determine_threads
    shell:
        """
         python workflow/envs/figaro/figaro/figaro.py -i data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro -o data/favabean/{wildcards.batch}-{wildcards.region}/figaro/ -f 1 -r 1 -a {params.expected_length} -F illumina > {log} 2>&1
         # Do some clean up of files that were solely used for figaro parameters generation
         rm -rf data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt/filteredForFigaro
        """


# Note that this uses the cutadapt files, and not the figaro ones! Need to change it so that it can do that for each pair of sequences separately instead of per batch
# Use this as your input: touch("data/favabean/{batch}-{region}/.tmp/.{sample}-cutadapt.done")
rule dada2_1_filterTrim:
    input:
        "data/favabean/{batch}-{region}/figaro/trimParameters.json"
    output:
        trimFilter=directory("data/favabean/{batch}-{region}/dada2"),
        donefile  = touch("data/favabean/{batch}-{region}/.tmp/.DADA2_trimFilter.done")
    log:
        R1="data/logs/dada2-1trimFilter-{batch}-{region}.log"
    conda:
        "../envs/dada2.yaml"
    params:
        figaro=config["figaro"]
    message: "DADA2 - trimming and filtering FASTQ files based on Figaro's {params.figaro} chosen method. batch: {wildcards.batch}, region: {wildcards.region}"
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2_1filterAndTrim.R -i {input} -f data/favabean/{wildcards.batch}-{wildcards.region}/cutadapt -o data/favabean/{wildcards.batch}-{wildcards.region}/dada2 -p {params.figaro} -c {threads} >> {log} 2>&1
        """


rule dada2_2_learnErrors_R1:
    input:
        "data/favabean/{batch}-{region}/dada2"
    output:
        "data/favabean/{batch}-{region}/dada2_DADA2Errors-R1.RData"
    log:
        "data/logs/dada2-2learnErrors-{batch}-{region}-R1.log"
    conda:
        "../envs/dada2.yaml"
    params:
        R="R1"
    threads:
        all_threads
    message: "DADA2 - Learning errors for R1. batch: {wildcards.batch}, region: {wildcards.region}"
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2_2LearnErrors.R -i {input} -r {params.R} -c {threads} >> {log} 2>&1
        """

use rule dada2_2_learnErrors_R1 as dada2_2_learnErrors_R2 with:
    input:
        "data/favabean/{batch}-{region}/dada2"
    output:
        "data/favabean/{batch}-{region}/dada2_DADA2Errors-R2.RData"
    log:
        "data/logs/dada2-2learnErrors-{batch}-{region}-R2.log"
    conda:
        "../envs/dada2.yaml"
    params:
        R="R2"
    message: "DADA2 - Learning errors for R2. batch: {wildcards.batch}, region: {wildcards.region}"



rule dada2_3_denoise_R1:
    input:
        "data/favabean/{batch}-{region}/dada2_DADA2Errors-R1.RData"
    output:
        "data/favabean/{batch}-{region}/dada2_DADA2Denoise-R1.RData"
    log:
        "data/logs/dada2-3denoise-{batch}-{region}-R1.log"
    conda:
        "../envs/dada2.yaml"
    params:
        R="R1"
    message: "DADA2 - Denoising amplicons based on the learned errors for R1 to create ASVs. batch: {wildcards.batch}, region: {wildcards.region}"
    threads:
        all_threads
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2_3denoise.R -i {input} -c {threads} -r {params.R} >> {log} 2>&1
        """

use rule dada2_3_denoise_R1 as dada2_3_denoise_R2 with:
    input:
        "data/favabean/{batch}-{region}/dada2_DADA2Errors-R2.RData"
    output:
        "data/favabean/{batch}-{region}/dada2_DADA2Denoise-R2.RData"
    log:
        "data/logs/dada2-3denoise-{batch}-{region}-R2.log"
    conda:
        "../envs/dada2.yaml"
    threads:
        determine_threads
    params:
        R="R2"
    message: "DADA2 - Denoising amplicons based on the learned errors for R2 to create ASVs. batch: {wildcards.batch}, region: {wildcards.region}"

rule dada2_4_mergePairedEnds:
    input:
        R1="data/favabean/{batch}-{region}/dada2_DADA2Denoise-R1.RData",
        R2="data/favabean/{batch}-{region}/dada2_DADA2Denoise-R2.RData"
    output:
        "data/favabean/{batch}-{region}/seqtab.tsv"
    log:
        "data/logs/dada2-mergePairedEnds-{batch}-{region}.log"
    conda:
        "../envs/dada2.yaml"
    message: "DADA2 - Merging ASV paired ends (R1, R2). batch: {wildcards.batch}, region: {wildcards.region}"
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2_4mergeR1R2.R -i {input.R1} -s {input.R2} >> {log} 2>&1
         rm -rf data/favabean/{wildcards.batch}-{wildcards.region}/.tmp
        """

combinations = [(row['Batch_ID'], row['region']) for _, row in samples_table.drop_duplicates(['Batch_ID', 'region']).iterrows()]

rule all:
    input:
        expand(
            [
                "data/favabean/{batch}-{region}/figaro/trimParameters.json",
                "data/favabean/{batch}-{region}/dada2",
                "data/favabean/{batch}-{region}/.tmp/.DADA2_trimFilter.done",
                # Include other target outputs as needed
            ],
            zip,
            batch=[x[0] for x in combinations],
            region=[x[1] for x in combinations]
        )


rule dada2_5_ChimeraDetectAndRemove:
    input:
        expand(
            "data/favabean/{batch}-{region}/seqtab.tsv",
            zip,
            batch=[combo[0] for combo in combinations],
            region=[combo[1] for combo in combinations]
        )
    output:
        "data/favabean/{region}_chimeraRemoved.RObjects"
    log:
        "data/logs/dada2-ChimeraASVs-{region}.log"
    conda:
        "../envs/dada2.yaml"
    threads:
        all_threads
    message: "DADA2 - Combining the results from the different runs, then removing chimeras through a de novo process. region: {wildcards.region}"
    shell:
        """
        Rscript --vanilla workflow/scripts/dada2_5chimera.R -i data/favabean/ -p {wildcards.region} -c {threads} -o {output} >> {log} 2>&1
        """


rule dada2_6_condense:
    input:
        "data/favabean/{region}_chimeraRemoved.RObjects"
    output:
        "data/favabean/{region}_condense.tsv"
    log:
        "data/logs/dada2-condenseASVs-{region}.log"
    conda:
        "../envs/dada2.yaml"
    threads:
        determine_threads
    message: "DADA2 - Condensing ASVs. region: {wildcards.region}"
    shell:
        """
         Rscript --vanilla workflow/scripts/dada2_6condense.R -i {input} -o {output} -r {wildcards.region} >> {log} 2>&1
        """

def get_urls(db_name, filetype):
    if config["taxonomy_database"][db_name]["use"]:
        if filetype == "ref":
            return config["taxonomy_database"][db_name]["url"]
        elif filetype == "species":
            return config["taxonomy_database"][db_name]["species"]
    return None

rule download_taxonomy_databases:
    output:
        ref     = "data/resources/{db}_ref.fa.gz",
        species = "data/resources/{db}_species.fa.gz"
    params:
        ref_url = lambda wildcards: get_urls(wildcards.db, "ref"),
        species_url = lambda wildcards: get_urls(wildcards.db, "species")
    run:
        if params.ref_url:
            shell("curl -L {params.ref_url} -o {output.ref}")
        if params.species_url:
            shell("curl -L {params.species_url} -o {output.species}")

rule dada2_7_assignTaxonomy:
    input:
        ASVs="data/favabean/{region}_condense.tsv",
        ref =ancient("data/resources/{db}_ref.fa.gz"),
        spec=ancient("data/resources/{db}_species.fa.gz")
    output:
        OTU_table="data/favabean/{region}_{db}_OTU.tsv",
        taxonomy ="data/favabean/{region}_{db}_taxonomy.tsv"
    log:
        "data/logs/dada2-{region}-{db}_taxonomy.log"
    conda:
        "../envs/dada2.yaml"
    threads:
        determine_threads
    message: "DADA2 - Assigning taxonomy to ASVs. region: {wildcards.region} using {wildcards.db} database."
    shell:
        """
        Rscript --vanilla workflow/scripts/dada2_7assignTaxonomy.R \
            -i {input.ASVs} \
            -d {input.ref} \
            -s {input.spec} \
            -c {threads} \
            -o {output.OTU_table} \
            -t {output.taxonomy} > {log} 2>&1
        """

rule paired_taxonomy:
    input:
        expand("data/favabean/{region}_{db}_OTU.tsv",region=[combo[1] for combo in combinations],db=[db for db in config["taxonomy_database"] if config["taxonomy_database"][db].get("use", False)])
    output:
        touch("data/favabean/.doneprimer_averaged.txt")
    log:
        "data/logs/primer_averaging.log"
    conda:
        "../envs/biom.yaml"
    message: "Creating biom files, and primer averaging, if needed."
    shell:
        """
        Rscript --vanilla workflow/scripts/primer_average.R > {log} 2>&1
        ls data/favabean/*OTU.tsv | parallel 'biom convert -i {{}} -o {{.}}.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy' > {log} 2>&1
        ls data/favabean/*taxonomy.tsv | parallel 'biom convert -i {{}} -o {{.}}.biom --to-json --table-type="OTU table"' > {log} 2>&1
        if [ -f data/favabean/primer_averaged.tsv ]; then
          ls data/favabean/primer_averaged.tsv | parallel 'biom convert -i {{}} -o {{.}}.biom --to-json --table-type="OTU table"' > {log} 2>&1
        else
          echo "File data/favabean/primer_averaged.tsv does not exist, skipping step."
        fi
        """

rule paired:
    input:
        expand("data/favabean/{region}_condense.tsv",region=[combo[1] for combo in combinations])
