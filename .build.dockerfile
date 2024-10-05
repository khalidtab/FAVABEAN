FROM snakemake/snakemake:v7.32.4
WORKDIR ../../
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache

RUN git clone https://github.com/khalidtab/FAVABEAN && mv ./FAVABEAN/workflow/build.snakefile ./build.snakefile && cp -r FAVABEAN/workflow ./workflow && cp FAVABEAN/workflow/snakefile ./snakefile && cp FAVABEAN/files_info.csv ./files_info.csv && cp FAVABEAN/favabean.yaml ./favabean.yaml && rm -r FAVABEAN && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile ./build.snakefile 2>&1 | tee -a environments.txt && rm ./build.snakefile && conda clean -a -y && rm ./snakefile

COPY start.sh /usr/local/bin/start
RUN chmod +x /usr/local/bin/start
CMD ["/usr/local/bin/start"]
