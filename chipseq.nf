#!/usr/bin/env nextflow

include { TRIM_GALORE } from './modules/common/trim_galore.nf'
include { BOWTIE2_ALIGN } from './modules/common/bowtie2_align.nf'
include { DEDUPLICATE } from './modules/common/deduplicate.nf'
include { COUNT_FRAGMENTS } from './modules/common/count_fragments.nf'
include { AVERAGE_REPLICATES } from './modules/cutandrun/average_replicates.nf'
include { SEACR_CALLPEAKS } from './modules/cutandrun/seacr_callpeaks.nf'
include { BIGWIG_FOLD } from './modules/cutandrun/bigwig_fold.nf'
include { MACS2_CALLPEAKS } from './modules/cutandrun/macs2_callpeaks_control.nf'
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
            [[index:++index, sample:row.sample, control:row.control, genome:row.genome, replicate:row.replicate], file("data/"+row.fastq_1), file("data/"+row.fastq_2)]
        }
    
    // Main pipeline
    TRIM_GALORE(samplesheet_ch)
    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed_reads)
    DEDUPLICATE(BOWTIE2_ALIGN.out.bam)
    COUNT_FRAGMENTS(DEDUPLICATE.out.nodupbam)

    grouped_bedgraph_ch = COUNT_FRAGMENTS.out.bedgraph
        .map { meta, bedgraph ->
            tuple (meta.subMap(['sample', 'control', 'genome']), bedgraph)
        }.groupTuple()
    
    AVERAGE_REPLICATES(grouped_bedgraph_ch)

    matched_bedgraph_ch = AVERAGE_REPLICATES.out.bedgraph
            .map { meta, bedgraph -> [meta.control, meta, bedgraph]}
        .combine(AVERAGE_REPLICATES.out.bedgraph
            .map { meta, bedgraph -> [meta.sample, meta, bedgraph]}, by:0)
        .map { key, meta1, bedgraph1, meta2, bedgraph2 ->
            [samplemeta:meta1, samplebedgraph:bedgraph1, controlbedgraph:bedgraph2]
        }

    SEACR_CALLPEAKS(matched_bedgraph_ch)

    matched_bigwig_ch = AVERAGE_REPLICATES.out.bigwig
            .map { meta, bigwig -> [meta.control, meta, bigwig]}
        .combine(AVERAGE_REPLICATES.out.bigwig
            .map { meta, bigwig -> [meta.sample, meta, bigwig]}, by:0)
        .map { key, meta1, bigwig1, meta2, bigwig2 ->
            [samplemeta:meta1, samplebigwig:bigwig1, controlbigwig:bigwig2]
        }

    BIGWIG_FOLD(matched_bigwig_ch)

    matched_nodupbam_ch = DEDUPLICATE.out.nodupbam
            .map { meta, nodupbam -> [meta.control, meta, nodupbam]}
        .combine(DEDUPLICATE.out.nodupbam
            .map { meta, nodupbam -> [meta.sample, meta, nodupbam]}, by:0)
        .filter(it -> it[3].replicate == "1")
        .map { key, meta1, nodupbam1, meta2, nodupbam2 ->
            [samplemeta:meta1, samplebam:nodupbam1, controlbam:nodupbam2]
        }
    
    MACS2_CALLPEAKS(matched_nodupbam_ch)

    SEACR_CALLPEAKS.out.seacrbed
        .concat(MACS2_CALLPEAKS.out.macs2bed)
        .filter { meta, bed -> bed.size() > 0}
        .set { bed }

    BED_TO_BIGBED(bed)

    // Compile QC summary and UCSC tracks
    matched_logs_ch = BOWTIE2_ALIGN.out.bowtie2log
        .combine(DEDUPLICATE.out.deduplogs, by:0)
    
    EXTRACT_QC(matched_logs_ch)

    EXTRACT_QC.out.summary.collectFile(name: 'qc_summary.txt', newLine: false, storeDir: "results", sort: { row -> row.splitCsv(sep:"\t")[0][0].toInteger() })
    COUNT_FRAGMENTS.out.bigwigtrack.collectFile(name: 'bigwigtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })
    AVERAGE_REPLICATES.out.avgbigwigtrack.collectFile(name: 'avgbigwigtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })
    BIGWIG_FOLD.out.fcbigwigtrack.collectFile(name: 'fcbigwigtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })
    BED_TO_BIGBED.out.bigbedtrack.collectFile(name: 'bigbedtracks.txt', newLine: false, storeDir: "results/ucsctracks", sort: { row -> row.splitCsv(sep:" ")[0][2] })

}
