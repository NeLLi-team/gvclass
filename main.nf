#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Pipeline version
params.version = '1.0.0'

// Default parameters
params.querydir = 'example'
params.outdir = "${params.querydir}_results"
params.database_path = "${workflow.projectDir}/resources"
params.keep_temp = true
params.mafftoption = 'linsi'
params.treeoption = 'iqtree'
params.mode_fast = true

// Log information
log.info """
GVClass Pipeline v${params.version}
===================================
querydir       : ${params.querydir}
outdir         : ${params.outdir}
database_path  : ${params.database_path}
keep_temp      : ${params.keep_temp}
mafft option   : ${params.mafftoption}
tree option    : ${params.treeoption}
mode_fast      : ${params.mode_fast}
"""

// Include modules
include { LOG_INFO } from './modules/log_info'
include { REFORMAT_FAA } from './modules/reformat'
include { GENECALLING } from './modules/genecalling'
include { RUN_HMMSEARCH } from './modules/hmmsearch'
include { EXTRACT_QHITS } from './modules/extract_qhits'
include { BLASTP_REDUCE_MERGE } from './modules/blastp'
include { ALIGN_TRIM } from './modules/align_trim'
include { BUILD_TREES } from './modules/build_trees'
include { GET_NN } from './modules/get_nn'
include { SUMMARIZE } from './modules/summarize'
include { COMBINEDOUT } from './modules/combinedout'
include { CLEANUP } from './modules/cleanup'

// Check and download resources if needed
process CHECK_RESOURCES {
    output:
    path 'resources_checked.txt', emit: checked

    script:
    """
    python ${workflow.projectDir}/bin/check_resources.py \
        --database_path ${params.database_path} \
        --output resources_checked.txt
    """
}

// Main workflow
workflow {
    // Check resources
    CHECK_RESOURCES()

    // Log run information
    LOG_INFO()

    // Find input files
    fna_files = Channel.fromPath("${params.querydir}/*.fna", checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    faa_files = Channel.fromPath("${params.querydir}/*.faa", checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    // Process FNA files (gene calling)
    GENECALLING(fna_files, CHECK_RESOURCES.out.checked)

    // Process FAA files (reformatting)
    REFORMAT_FAA(faa_files, CHECK_RESOURCES.out.checked)

    // Combine outputs from gene calling and reformatting
    query_faa = GENECALLING.out.faa.mix(REFORMAT_FAA.out.faa)
    query_stats = GENECALLING.out.stats.mix(REFORMAT_FAA.out.stats)

    // Run HMMSearch on all query FAA files
    RUN_HMMSEARCH(query_faa, CHECK_RESOURCES.out.checked)

    // Get model names from HMM file
    process GET_MODEL_NAMES {
        input:
        path resources_checked

        output:
        path 'model_names.txt', emit: model_names

        script:
        """
        python ${workflow.projectDir}/bin/get_model_names.py \
            --hmm_file ${params.database_path}/models/combined.hmm \
            --mode_fast ${params.mode_fast} \
            --output model_names.txt
        """
    }

    // Execute the GET_MODEL_NAMES process
    GET_MODEL_NAMES(CHECK_RESOURCES.out.checked)

    // Extract query hits for each model
    model_names = GET_MODEL_NAMES.out.model_names
        .splitText()
        .map { it.trim() }

    // Create combinations of query and model
    query_model_combinations = RUN_HMMSEARCH.out.filtered
        .combine(model_names)
        .map { querybase, filtered, query_faa, model -> tuple(querybase, filtered, query_faa, model) }

    // Extract hits for each model
    EXTRACT_QHITS(query_model_combinations)

    // BLASTP, alignment, and tree building for each query-model combination
    BLASTP_REDUCE_MERGE(EXTRACT_QHITS.out)
    ALIGN_TRIM(BLASTP_REDUCE_MERGE.out.merged)
    BUILD_TREES(ALIGN_TRIM.out.trimmed)

    // Group tree files by query
    tree_files = BUILD_TREES.out.tree
        .map { querybase, modelbase, tree -> tuple(querybase, tree) }
        .groupTuple()

    // Group alignment files by query
    aln_files = ALIGN_TRIM.out.trimmed
        .map { querybase, modelbase, aln -> tuple(querybase, aln) }
        .groupTuple()

    // Get nearest neighbors
    GET_NN(tree_files.join(aln_files))

    // Summarize results
    SUMMARIZE(
        GET_NN.out.nn
            .join(RUN_HMMSEARCH.out.counts)
            .join(query_stats)
    )

    // Combine all summaries
    all_summaries = SUMMARIZE.out.summary.map { querybase, summary -> summary }.collect()
    COMBINEDOUT(all_summaries)

    // Force a clean run of the COMBINEDOUT process
    process FORCE_COMBINEDOUT {
        publishDir "${params.outdir}", mode: 'copy'

        input:
        path summaries

        output:
        path "gvclass_out_v${params.version}.tab", emit: combined

        script:
        """
        # Create a directory structure for each summary file
        for summary in ${summaries.join(' ')}; do
            filename=\$(basename \$summary)
            dirname=\$(echo \$filename | cut -d '.' -f 1)
            mkdir -p temp_summaries/\$dirname
            cp \$summary temp_summaries/\$dirname/
        done

        python ${workflow.projectDir}/bin/combinedout.py \\
            -r temp_summaries \\
            -o gvclass_out_v${params.version}.tab \\
            -v ${params.version}
        """
    }

    FORCE_COMBINEDOUT(all_summaries)

    // Cleanup and compress results
    SUMMARIZE.out.summary
        .combine(FORCE_COMBINEDOUT.out.combined)
        .map { querybase, summary, combined -> tuple(querybase, summary, combined) }
        .set { cleanup_input }

    CLEANUP(cleanup_input)
}
