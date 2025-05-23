process CLEANUP {
    tag "${querybase}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(querybase), path(summary), path(combined)

    output:
    path "${querybase}.tar.gz"

    script:
    """
    # Create a temporary directory structure
    mkdir -p temp/${querybase}

    # Copy the summary file to the temporary directory
    cp ${summary} temp/${querybase}/

    # Create the tar.gz file
    tar -zcvf ${querybase}.tar.gz -C temp ${querybase}

    # Clean up the temporary directory
    rm -rf temp
    """
}
