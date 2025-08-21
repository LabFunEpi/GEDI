#!/usr/bin/env nextflow

include { TRIM_GALORE } from './modules/common/trim_galore.nf'
include { BOWTIE2_ALIGN } from './modules/common/bowtie2_align.nf'
include { DEDUPLICATE } from './modules/common/deduplicate.nf'
include { COUNT_FRAGMENTS } from './modules/common/count_fragments.nf'
include { MACS2_CALLPEAKS } from './modules/common/macs2_callpeaks.nf'
include { BED_TO_BIGBED } from './modules/common/bed_to_bigbed.nf'
include { EXTRACT_QC } from './modules/common/extract_qc.nf'

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
    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed_reads)
    DEDUPLICATE(BOWTIE2_ALIGN.out.bam)
    COUNT_FRAGMENTS(DEDUPLICATE.out.nodupbam)
    MACS2_CALLPEAKS(DEDUPLICATE.out.nodupbam)
    BED_TO_BIGBED(MACS2_CALLPEAKS.out.macs2bed)

    // Compile QC summary and UCSC tracks
    matched_logs_ch = BOWTIE2_ALIGN.out.bowtie2log
        .combine(DEDUPLICATE.out.deduplogs, by:0)
    
    EXTRACT_QC(matched_logs_ch)

    EXTRACT_QC.out.summary.collectFile(name: 'qc_summary.txt', newLine: false, storeDir: "results", sort: { row -> row.splitCsv(sep:"\t")[0][0].toInteger() })
    COUNT_FRAGMENTS.out.bigwigtrack.collectFile(name: 'bigwigtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })
    BED_TO_BIGBED.out.bigbedtrack.collectFile(name: 'bigbedtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })

}
