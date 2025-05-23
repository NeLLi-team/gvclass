process COMBINEDOUT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path summaries

    output:
    path "gvclass_out_v${params.version}.tab", emit: combined

    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p ${params.outdir}

    # Create a directory structure for each summary file
    for summary in ${summaries.join(' ')}; do
        filename=\$(basename \$summary)
        dirname=\$(echo \$filename | cut -d '.' -f 1)
        mkdir -p temp_summaries/\$dirname
        cp \$summary temp_summaries/\$dirname/
    done

    python ${workflow.projectDir}/bin/combinedout.py \
        -r temp_summaries \
        -o gvclass_out_v${params.version}.tab \
        -v ${params.version}
    """
}
