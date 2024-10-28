QUERYDIR=$1
PROCESSES=$2

apptainer run --containall --bind $(pwd):/workdir --pwd /workdir  \
  docker://docker.io/doejgi/gvclass:latest \
  snakemake --snakefile /gvclass/workflow/Snakefile \
           -j $PROCESSES \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir="/workdir/$QUERYDIR" \
           database_path="/gvclass/resources"
