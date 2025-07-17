#!/usr/bin/env nextflow

process AVERAGE_REPLICATES {

    publishDir "results/bedgraphs", mode: 'copy', pattern: "${meta.sample}.bg"
    publishDir "results/bigwigs", mode: 'copy', pattern: "${meta.sample}.bw"

    input:
    tuple val(meta), path(bedgraphs)

    output:
    tuple val(meta), path("${meta.sample}.bg"), emit: bedgraph
    tuple val(meta), path("${meta.sample}.bw"), emit: bigwig
    path("${meta.sample}.bigwigtrack.txt"), emit: avgbigwigtrack

    script:
    """
    
    module load libbigwig
    module load htslib
    module load gsl
    module load wiggletools
    module load ucscbigtools
    
    wiggletools write_bg ${meta.sample}.bg mean ${bedgraphs}
    
    bedGraphToBigWig ${meta.sample}.bg ${projectDir}/references/chromsizes/${meta.genome} ${meta.sample}.bw
    
    printf 'track type=bigWig name=\"${meta.sample}_bw\" bigDataUrl=cloudpath/${meta.sample}.bw description=${meta.sample} visibility=full smoothingWindow=4 windowingFunction=mean maxHeightPixels=100:50:8\\n' > ${meta.sample}.bigwigtrack.txt
    
    """
}
