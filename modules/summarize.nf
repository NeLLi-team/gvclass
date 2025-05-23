process SUMMARIZE {
    tag "${querybase}"
    publishDir "${params.outdir}/${querybase}", mode: 'copy'

    input:
    tuple val(querybase), path(nn_tree), path(models_count), path(querystats)

    output:
    tuple val(querybase), path("${querybase}.summary.tab"), emit: summary

    script:
    """
    python ${workflow.projectDir}/bin/summarize.py \
        -n ${nn_tree} \
        -g ${models_count} \
        -q ${querystats} \
        -s ${querybase}.summary.tab \
        -f ${params.database_path}/order_completeness.tab
    """
}
