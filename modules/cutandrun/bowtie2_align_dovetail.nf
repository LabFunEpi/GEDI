#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {

    // publishDir "results", mode: 'copy', pattern: "${prefix}.bam*"
    publishDir "results/qc/bowtie2_reports", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}.bowtie2.log"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sample}_${meta.replicate}.bam"), emit: bam
    path "${meta.sample}_${meta.replicate}.bam.bai", emit: bai
    tuple val(meta), path("${meta.sample}_${meta.replicate}.bowtie2.log"), emit: bowtie2log

    script:
    """
    
    module load bowtie2
    module load samtools
    
    bowtie2 --dovetail --threads $task.cpus -x ${projectDir}/references/bowtie2/${meta.genome} -1 ${fastq_1} -2 ${fastq_2} 2> ${meta.sample}_${meta.replicate}.bowtie2.log | \
        samtools view -b -q 30 - | samtools sort -O BAM - > ${meta.sample}_${meta.replicate}_temp.bam
    samtools index ${meta.sample}_${meta.replicate}_temp.bam
    samtools view -b ${meta.sample}_${meta.replicate}_temp.bam `seq 1 22 | sed 's/^/chr/'` chrX chrY > ${meta.sample}_${meta.replicate}.bam
    samtools index ${meta.sample}_${meta.replicate}.bam
    
    """
}
