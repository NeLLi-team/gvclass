QUERYDIR=$1
PROCESSES=$2

docker run \
  -v $(pwd):$(pwd) -w $(pwd) doejgi/gvclass:latest \
  snakemake --snakefile /gvclass/workflow/Snakefile \
           -j $PROCESSES \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir="$QUERYDIR" \
           database_path="/gvclass/resources"
