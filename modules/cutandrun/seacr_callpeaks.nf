#!/usr/bin/env nextflow

process SEACR_CALLPEAKS {

    publishDir "results/beds", mode: 'copy', pattern: "${meta.sample}.seacr.bed"

    input:
    tuple val(meta), path(samplebedgraph), path(controlbedgraph)

    output:
    tuple val(meta), path("${meta.sample}.seacr.bed"), emit: seacrbed

    script:
    """
    
    module load r
    module load bedtools
    
    /research/bsi/tools/biotools/seacr/1.3/SEACR-1.3/SEACR_1.3.sh ${samplebedgraph} ${controlbedgraph} norm relaxed ${meta.sample}

    mv "${meta.sample}.relaxed.bed" "${meta.sample}.seacr.bed"
    
    """
}
