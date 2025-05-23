process GENECALLING {
    tag "${querybase}"
    publishDir "${params.outdir}/${querybase}", mode: 'copy', pattern: "*.{gff,tab,faa,fna}"

    input:
    tuple val(querybase), path(fna_file)
    path resources_checked

    output:
    tuple val(querybase), path("query_faa/${querybase}.faa"), emit: faa
    tuple val(querybase), path("stats/${querybase}.stats.tab"), emit: stats
    path "query_gff/${querybase}.gff"
    path "stats/${querybase}.genecalling.tab"
    path "query_fna/${querybase}.fna"
    path "hmmout/models.out", optional: true

    script:
    """
    mkdir -p query_gff stats query_faa query_fna hmmout

    if [ -s ${fna_file} ]; then
        python ${workflow.projectDir}/bin/opgecall_no_click.py \
            -f ${fna_file} \
            -g query_gff/${querybase}.gff \
            -gs stats/${querybase}.genecalling.tab \
            -ss stats/${querybase}.stats.tab \
            -fa query_faa/${querybase}.faa \
            -m ${params.database_path}/models/combined.hmm \
            -fn query_fna/${querybase}.fna \
            -o hmmout/models.out
    else
        echo "Input file ${fna_file} is not available. Skipping gene calling."
        touch query_gff/${querybase}.gff stats/${querybase}.genecalling.tab \
              stats/${querybase}.stats.tab query_fna/${querybase}.fna \
              hmmout/models.out query_faa/${querybase}.faa
    fi
    """
}
