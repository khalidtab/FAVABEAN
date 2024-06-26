FROM snakemake/snakemake:v7.32.4
WORKDIR ../../
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/khalidtab/FAVABEAN && mv ./FAVABEAN/workflow/build.snakefile ./build.snakefile && cp -r FAVABEAN/workflow ./workflow && cp FAVABEAN/workflow/snakefile ./snakefile && cp FAVABEAN/files_info.csv ./files_info.csv && cp FAVABEAN/input.yaml ./input.yaml && rm -r FAVABEAN && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile ./build.snakefile > environments.txt && rm ./build.snakefile && conda clean -a -y && rm ./snakefile
