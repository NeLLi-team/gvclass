process ALIGN_TRIM {
    tag "${querybase}:${modelbase}"
    publishDir "${params.outdir}/${querybase}/queryrefs_aligned", mode: 'copy', enabled: params.keep_temp

    input:
    tuple val(querybase), val(modelbase), path(mergedfaa)

    output:
    tuple val(querybase), val(modelbase), path("${modelbase}.mafft"), emit: aligned
    tuple val(querybase), val(modelbase), path("${modelbase}.mafft01"), emit: trimmed

    script:
    """
    python ${workflow.projectDir}/bin/align_trim.py \
        -q ${mergedfaa} \
        -a ${modelbase}.mafft \
        -t ${modelbase}.mafft01 \
        -o ${params.mafftoption}

    touch ${modelbase}.mafft ${modelbase}.mafft01
    """
}
