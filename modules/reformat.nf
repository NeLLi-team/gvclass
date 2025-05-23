process REFORMAT_FAA {
    tag "${querybase}"
    publishDir "${params.outdir}/${querybase}", mode: 'copy', pattern: "*.{faa,tab}"

    input:
    tuple val(querybase), path(faa_file)
    path resources_checked

    output:
    tuple val(querybase), path("query_faa/${querybase}.faa"), emit: faa
    tuple val(querybase), path("stats/${querybase}.stats.tab"), emit: stats

    script:
    """
    mkdir -p query_faa stats

    python ${workflow.projectDir}/bin/reformat_no_click.py \
        -i ${faa_file} \
        -o query_faa/${querybase}.faa \
        -s stats/${querybase}.stats.tab
    """
}
