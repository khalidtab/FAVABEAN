FROM snakemake/snakemake:v7.32.4
WORKDIR ../../
RUN curl -s https://www.random.org/integers/?num=1&min=1&max=100000&col=1&base=10&format=plain&rnd=new > random_number.txt
RUN git clone https://github.com/khalidtab/FAVABEAN && mv ./FAVABEAN/workflow/build.snakefile ./build.snakefile && cp -r FAVABEAN/workflow ./workflow && cp FAVABEAN/workflow/snakefile ./snakefile && rm -r FAVABEAN && snakemake snakefile.final --use-conda --cores all --conda-create-envs-only --snakefile ./build.snakefile > environments.txt && rm ./build.snakefile && conda clean -a -y && rm ./snakefile

