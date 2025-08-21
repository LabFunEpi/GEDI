#!/usr/bin/env nextflow

process STAR_ALIGN {

    // publishDir "results", mode: 'copy', pattern: "${prefix}.bam*"
    publishDir "results/qc/star_reports", mode: 'copy', pattern: "*.out"
    publishDir "results/counts", mode: 'copy', pattern: "*.tab"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sample}_${meta.replicate}.bam"), emit: bam
    path "${meta.sample}_${meta.replicate}.bam.bai", emit: bai
    path("*.out")
    path("*.tab")

    script:
    """
    
    module load star
    module load samtools

    STAR --runMode alignReads \\
        --runThreadN $task.cpus \\
        --genomeDir ${params.refDir}/star/${meta.genome} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --outFileNamePrefix ${meta.sample}_${meta.replicate}_ \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNmax 10 \\
        --outFilterType BySJout \\
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --readFilesCommand zcat
    
    samtools index ${meta.sample}_${meta.replicate}_Aligned.sortedByCoord.out.bam

    samtools view -b ${meta.sample}_${meta.replicate}_Aligned.sortedByCoord.out.bam \\
            `seq 1 22 | sed 's/^/chr/'` chrX chrY | \\
        samtools view -b -q 30 - | \\
        samtools sort -O BAM - > ${meta.sample}_${meta.replicate}.bam
    
    samtools index ${meta.sample}_${meta.replicate}.bam
    
    """
}
