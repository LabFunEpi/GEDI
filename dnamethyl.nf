#!/usr/bin/env nextflow

include { TRIM_GALORE } from './modules/common/trim_galore.nf'
// include { BISMARK_ALIGN } from './modules/dnamethyl/bismark_align.nf'
// include { BISMARK_DEDUP } from './modules/dnamethyl/bismark_dedup.nf'
// include { BISMARK_EXTRACT } from './modules/dnamethyl/bismark_extract.nf'

// Primary input
params.input_csv = "data/samples.csv"

workflow {

    // Read and parse input
    def index = 0
    samplesheet_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row ->
            [[index:++index, sample:row.sample, genome:row.genome, replicate:row.replicate], file("data/"+row.fastq_1), file("data/"+row.fastq_2)]
        }
    
    // Main pipeline
    TRIM_GALORE(samplesheet_ch)
    // BISMARK_ALIGN(TRIM_GALORE.out.trimmed_reads)
    // BISMARK_DEDUP(BISMARK_ALIGN.out.bam)
    // BISMARK_EXTRACT(BISMARK_DEDUP.out.nodupbam)

}
