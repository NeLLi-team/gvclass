QUERYDIR=$1
PROCESSES=$2

shifterimg pull docker:doejgi/gvclass:latest

shifter --image=docker:doejgi/gvclass:latest  \
  snakemake --snakefile /gvclass/workflow/Snakefile \
           -j $PROCESSES \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir="$QUERYDIR" \
           database_path="/gvclass/resources"
