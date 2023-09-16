rule alpha_div_calc: # Provides per sample alpha calculation
   version: "1.0"
   conda:
      "../../workflow/envs/cutadapt.yaml"
   input:
      "data/tsv/{sample}.tsv"
   output:
      "data/alpha_div/calc_{sample}–{alpha}.txt"
   message: "Cutadapt for files"
   shell:
      "cutadapt -g {foward} -G {reverse} -o {forward.output} -p {reverse.output} {forward.input} {reverse.input} --minimum-length {minimum.length} >> {log}"

rule alpha:
  input:
   expand("tmp/expand–{sample}–{alpha}–{group}.txt", sample=config["mysample"], alpha=config["alpha"], group=config["group"])