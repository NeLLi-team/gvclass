process GET_NN {
    tag "${querybase}"
    publishDir "${params.outdir}/${querybase}/stats", mode: 'copy'

    input:
    tuple val(querybase), path(tree_files), path(aln_files)

    output:
    tuple val(querybase), path("${querybase}.tree_nn"), emit: nn

    script:
    """
    # Create directories to match expected structure
    mkdir -p queryrefs_genetrees queryrefs_aligned

    # Link tree files to expected location
    for tree_file in ${tree_files}; do
        ln -sf \$PWD/\$tree_file queryrefs_genetrees/
    done

    # Link alignment files to expected location
    for aln_file in ${aln_files}; do
        ln -sf \$PWD/\$aln_file queryrefs_aligned/
    done

    python ${workflow.projectDir}/bin/get_nn_tree.py \
        -q ${querybase} \
        -t queryrefs_genetrees/ \
        -o ${querybase}.tree_nn \
        -l ${params.database_path}/ncldvApril24_labels.txt
    """
}
