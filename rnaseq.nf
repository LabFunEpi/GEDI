#!/usr/bin/env nextflow

include { TRIM_GALORE } from './modules/common/trim_galore.nf'
include { STAR_ALIGN } from './modules/rnaseq/star_align.nf'
include { DEDUPLICATE } from './modules/common/deduplicate.nf'
include { HTSEQ_COUNT } from './modules/rnaseq/htseq_count.nf'

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
    STAR_ALIGN(TRIM_GALORE.out.trimmed_reads)
    DEDUPLICATE(STAR_ALIGN.out.bam)
    HTSEQ_COUNT(DEDUPLICATE.out.nodupbam)

}
