#!/usr/bin/env nextflow

process EXTRACT_QC {

    input:
    tuple val(meta), path(bowtie2log), path(picardlog), path(lclog)

    output:
    path("${meta.sample}_${meta.replicate}.qc.txt"), emit: summary

    script:
    """
    
    paste -d '\\t' \
        <(echo ${meta.index}) \
        <(echo ${meta.sample}_${meta.replicate}) \
        <(tail -1 ${bowtie2log} | cut -d"%" -f1) \
        <(awk -v OFS='\t' '{print \$1, \$5, \$6, \$7}' ${lclog}) \
        <(sed -n '8p' ${picardlog} | cut -f9) \
    > ${meta.sample}_${meta.replicate}.qc.txt
    
    """
}
