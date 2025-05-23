process EXTRACT_QHITS {
    tag "${querybase}:${modelbase}"
    publishDir "${params.outdir}/${querybase}/query_hits_faa", mode: 'copy', enabled: params.keep_temp

    input:
    tuple val(querybase), path(hmmout), path(queryfaa), val(modelbase)

    output:
    tuple val(querybase), val(modelbase), path("${modelbase}.faa")

    script:
    """
    if [ -s ${hmmout} ]; then
        python ${workflow.projectDir}/bin/extract_qhits.py \
            -h ${hmmout} \
            -q ${queryfaa} \
            -o ${modelbase}.faa
    else
        echo "No hits found in hmmsearch output. Skipping hit extraction."
        touch ${modelbase}.faa
    fi
    """
}
