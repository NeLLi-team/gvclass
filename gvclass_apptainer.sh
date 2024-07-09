QUERYDIR=$1

apptainer exec --writable-tmpfs \
  --bind ./$QUERYDIR:/gvclass/input \
  gvclass.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate snk && \
           cd /gvclass && \
           snakemake --snakefile workflow/Snakefile \
           -j 24 \
           --use-conda \
           --conda-frontend mamba \
           --conda-prefix /gvclass/.snakemake/conda \
           --config querydir='/gvclass/input' \
           database_path='/gvclass/resources'"
