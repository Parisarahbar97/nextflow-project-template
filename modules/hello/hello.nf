// Minimal example module (DSL2)

process HELLO_PROCESS {
    // Run this process inside a container image
    container 'ghcr.io/johnsonlab-ic/causal-flow:latest'

    tag { sampleName }
    publishDir "${params.outputDir}/${sampleName}", mode: 'copy'
    label 'low'

    input:
    val sampleName

    output:
    path "${sampleName}_greeting.txt"

    script:
    """
    echo "Hello from Nextflow, ${sampleName}!" > ${sampleName}_greeting.txt
    cat ${sampleName}_greeting.txt
    """
}

// NOTE: removed the workflow wrapper so this module only exposes the process
// Use the process directly from `main.nf` (e.g. include { HELLO_PROCESS } from './modules/hello/hello')