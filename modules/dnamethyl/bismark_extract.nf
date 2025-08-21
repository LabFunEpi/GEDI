#!/usr/bin/env nextflow

process HTSEQ_COUNT {

    publishDir "results/counts", mode: 'copy', pattern: "*.tab"

    input:
    tuple val(meta), path(bam)

    output:
    path("*.tab")

    script:
    """
    
    module load htseq
    
    htseq-count -f bam -s no -r pos ${bam} \\
        ${params.refDir}/gtf/${meta.genome}.gtf > \\
        ${meta.sample}_${meta.replicate}_htseq.tab
    
    """
}
