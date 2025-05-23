process LOG_INFO {
    publishDir "${params.outdir}", mode: 'copy'
    
    output:
    path 'run_info.log'
    
    script:
    """
    python ${workflow.projectDir}/bin/log_run_info.py \
        --version ${params.version} \
        --config '${groovy.json.JsonOutput.toJson(params)}' \
        --output run_info.log
    """
}
