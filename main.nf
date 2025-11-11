nextflow.enable.dsl=2

// include the hello process from the example module
include { HELLO_PROCESS } from './modules/hello/hello'

// Default parameters
params.samples = "samples.csv"
params.outputDir = "results"
params.help = false


workflow {
    if (params.help) {
        println "Usage: nextflow run main.nf --samples samples.csv"
        exit 0
    }

    // simple samples channel: read CSV with header `name` or fallback to built-in example
    def samples_ch = Channel
        .fromPath(params.samples)
        .ifEmpty { Channel.of([['name':'example']]) }
        .splitCsv(header:true)
        .map { row -> row.name }

    // call the HELLO_PROCESS for each sample name (process expects a channel of vals)
    hello_out = HELLO_PROCESS(samples_ch)

    // publish greetings — this will print the produced path(s)
    hello_out.view()
}