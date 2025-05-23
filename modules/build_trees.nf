process BUILD_TREES {
    tag "${querybase}:${modelbase}"
    publishDir "${params.outdir}/${querybase}/queryrefs_genetrees", mode: 'copy', enabled: params.keep_temp

    input:
    tuple val(querybase), val(modelbase), path(trimmedaln)

    output:
    tuple val(querybase), val(modelbase), path("${modelbase}.treefile"), emit: tree

    script:
    """
    python ${workflow.projectDir}/bin/build_tree.py \
        -a ${trimmedaln} \
        -t ${modelbase}.treefile \
        -m ${params.treeoption}

    touch ${modelbase}.treefile
    """
}
