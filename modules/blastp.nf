process BLASTP_REDUCE_MERGE {
    tag "${querybase}:${modelbase}"
    publishDir "${params.outdir}/${querybase}", mode: 'copy', enabled: params.keep_temp

    input:
    tuple val(querybase), val(modelbase), path(queryhitsfaa)

    output:
    tuple val(querybase), val(modelbase), path("blastp_out/${modelbase}.m8"), emit: blast
    tuple val(querybase), val(modelbase), path("query_hits_merged_faa/${modelbase}.faa"), emit: merged

    script:
    """
    mkdir -p blastp_out query_hits_merged_faa

    python ${workflow.projectDir}/bin/blastp_reduce_merge.py \
        -q ${queryhitsfaa} \
        -r ${params.database_path}/database/faa/${modelbase}.faa \
        -d ${params.database_path}/database/dmnd/${modelbase}.dmnd \
        -b blastp_out/${modelbase}.m8 \
        -o query_hits_merged_faa/${modelbase}.faa

    touch blastp_out/${modelbase}.m8 query_hits_merged_faa/${modelbase}.faa
    """
}
