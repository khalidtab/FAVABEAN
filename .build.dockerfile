FROM snakemake/snakemake:stable
WORKDIR ../../
RUN git clone https://github.com/khalidtab/FAVABEAN && mv ./FAVABEAN/workflow/build.snakefile ./build.snakefile && cp -r FAVABEAN/workflow ./workflow && cp FAVABEAN/workflow/snakefile ./snakefile && cp FAVABEAN/files_info.csv ./files_info.csv && cp FAVABEAN/input.yaml ./input.yaml && rm -r FAVABEAN && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile ./build.snakefile > environments.txt && rm ./build.snakefile && conda clean -a -y && rm ./snakefile
