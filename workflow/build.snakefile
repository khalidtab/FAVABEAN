# =============================================================================
# build.snakefile — Pre-create all conda environments at Docker build time
#
# Used by .build.dockerfile with:
#   snakemake snakefile.final --use-conda --cores all --conda-create-envs-only
#
# Each rule references one conda env yaml so that Snakemake resolves and
# caches every environment during the image build.
# =============================================================================

rule all:
    input:
        "cutadapt.done",
        "seqkit.done",
        "figaro.done",
        "dada2.done",
        "biom.done",
        "sidle.done"
    output:
        touch("snakefile.final")

rule build_cutadapt:
    output: touch("cutadapt.done")
    conda: "workflow/envs/cutadapt.yaml"
    shell: "echo 'cutadapt env created'"

rule build_seqkit:
    output: touch("seqkit.done")
    conda: "workflow/envs/seqkit.yaml"
    shell: "echo 'seqkit env created'"

rule build_figaro:
    output: touch("figaro.done")
    conda: "workflow/envs/figaro.yaml"
    shell: "echo 'figaro env created'"

rule build_dada2:
    output: touch("dada2.done")
    conda: "workflow/envs/dada2.yaml"
    shell: "echo 'dada2 env created'"

rule build_biom:
    output: touch("biom.done")
    conda: "workflow/envs/biom.yaml"
    shell: "echo 'biom env created'"

rule build_sidle:
    output: touch("sidle.done")
    conda: "workflow/envs/sidle.yaml"
    shell: "echo 'sidle env created'"
