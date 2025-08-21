#!/usr/bin/env nextflow

process BIGWIG_FOLD {

    publishDir "results/bigwigs", mode: 'copy', pattern: "${meta.sample}_FC.bw"

    input:
    tuple val(meta), path(samplebigwig), path(controlbigwig)

    output:
    tuple val(meta), path("${meta.sample}_FC.bw"), emit: fcbigwig
    path("${meta.sample}.bigwigtrack.txt"), emit: fcbigwigtrack

    script:
    """
    
    module load deeptools
    
    bigwigCompare --bigwig1 ${samplebigwig} --bigwig2 ${controlbigwig} --operation ratio --binSize 10 --outFileName ${meta.sample}_FC.bw --numberOfProcessors $task.cpus

    printf 'track type=bigWig name=\"${meta.sample}_FC_bw\" bigDataUrl=cloudpath/${meta.sample}_FC.bw description=${meta.sample} visibility=full smoothingWindow=4 windowingFunction=mean maxHeightPixels=100:50:8\\n' > ${meta.sample}.bigwigtrack.txt
    
    """
}
