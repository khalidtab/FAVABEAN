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

rule bbmap:
   conda:
      "./workflow/envs/bbmap.yaml"
   message: "Creating bbmap"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}bbmap"))
   shell:
      "touch {output}"

rule seqtk:
   conda:
      "./workflow/envs/seqtk.yaml"
   message: "Creating seqtk"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}seqtk"))
   shell:
      "touch {output}"


rule pigz:
   conda:
      "./workflow/envs/pigz.yaml"
   message: "Creating pigz"
   input:
      "{sample}"
   output:
      temporary(touch("{sample}pigz"))
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

rule results:
   version: "1.0"
   input:
      rules.figaro.output,
      rules.cutadapt.output,
      rules.bbmap.output,
      rules.seqtk.output,
      rules.pigz.output,
      rules.dada2.output
   message: "Cleaning up…"
   output:
      "{sample}.final"
   shell:
      "conda clean -a && echo Done building environments"