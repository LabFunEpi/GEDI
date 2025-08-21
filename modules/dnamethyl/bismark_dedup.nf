#!/usr/bin/env nextflow

process DEDUPLICATE {

    publishDir "results/bams", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}.NoDup.bam*"
    publishDir "results/qc/picard_reports", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}_picard.log"
    publishDir "results/qc/lc_reports", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}_LC.log"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_${meta.replicate}.NoDup.bam"), emit: nodupbam
    path "${meta.sample}_${meta.replicate}.NoDup.bam.bai", emit: nondupbai
    tuple val(meta), path("${meta.sample}_${meta.replicate}_picard.log"), path("${meta.sample}_${meta.replicate}_LC.log"), emit: deduplogs

    script:
    """
    
    module load picard
    module load samtools
    module load bedtools
    
    picard MarkDuplicates I=${bam} O=${meta.sample}_${meta.replicate}.NoDup.bam \\
        M=${meta.sample}_${meta.replicate}_picard.log REMOVE_DUPLICATES=true
    
    samtools index ${meta.sample}_${meta.replicate}.NoDup.bam
    
    picard MarkDuplicates I=${bam} O=DupMark.bam \\
        M=${meta.sample}_${meta.replicate}_picard_1.log REMOVE_DUPLICATES=false
    
    samtools index DupMark.bam
    
    bedtools bamtobed -i DupMark.bam | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,\$6}' | \\
        grep --color=auto -v 'chrM' | sort | uniq -c | \\
        awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0/mt,m1/m0,(m2==0)?0:m1/m2}' \\
        > ${meta.sample}_${meta.replicate}_LC.log
    
    """
}
