QUERYDIR=$1
PROCESSES=$2

apptainer run \
  docker://docker.io/doejgi/gvclass:latest \
  snakemake --snakefile /gvclass/workflow/Snakefile \
           -j $PROCESSES \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir="$QUERYDIR"
           database_path="/gvclass/resources"
