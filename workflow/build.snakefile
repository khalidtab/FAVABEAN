rule all:
   input: "{sample}.final"
   shell: "touch {input}"

rule figaro:
   conda:
      "./workflow/envs/figaro.yaml"
   message: "Creating Figaro"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}figaro"))
   shell:
      "touch {output}"

rule cutadapt:
   conda:
      "./workflow/envs/cutadapt.yaml"
   message: "Creating cutadapt"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}cutadapt"))
   shell:
      "touch {output}"

rule seqkit:
   conda:
      "./workflow/envs/seqkit.yaml"
   message: "Creating seqkit"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}seqkit"))
   shell:
      "touch {output}"


rule dada2:
   conda:
      "./workflow/envs/dada2.yaml"
   message: "Creating dada2"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}dada2"))
   shell:
      "touch {output}"

rule biom:
   conda:
      "./workflow/envs/biom.yaml"
   message: "Creating biom"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}biom"))
   shell:
      "touch {output}"


rule results:
   input:
      rules.figaro.output,
      rules.cutadapt.output,
      rules.seqkit.output,
      rules.dada2.output,
      rules.biom.output
   message: "Cleaning up…"
   output:
      "{sample}.final"
   shell:
      "conda clean -a && echo Done building environments"
