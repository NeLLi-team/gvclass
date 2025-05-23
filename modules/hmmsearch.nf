process RUN_HMMSEARCH {
    tag "${querybase}"
    publishDir "${params.outdir}/${querybase}/hmmout", mode: 'copy'

    input:
    tuple val(querybase), path(queryfaa)
    path resources_checked

    output:
    tuple val(querybase), path('models.out.filtered'), path(queryfaa), emit: filtered
    tuple val(querybase), path('models.counts'), emit: counts
    path 'models.score'
    path 'models.out'

    script:
    """
    if [ -s ${queryfaa} ]; then
        python ${workflow.projectDir}/bin/hmmsearch_no_click.py \
            -q ${queryfaa} \
            -m ${params.database_path}/models/combined.hmm \
            --hmmout models.out \
            -hf models.out.filtered \
            -c models.counts \
            -s models.score \
            -f ${params.database_path}/models_APRIL24--databaseApril24.cutoffs
    else
        echo "Input file ${queryfaa} is empty. Skipping hmmsearch."
        touch models.out models.counts models.score models.out.filtered
    fi
    """
}
