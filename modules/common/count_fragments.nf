#!/usr/bin/env nextflow

process COUNT_FRAGMENTS {

    publishDir "results/bedgraphs", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}.bg"
    publishDir "results/bigwigs", mode: 'copy', pattern: "${meta.sample}_${meta.replicate}.bw"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_${meta.replicate}.bg"), emit: bedgraph
    tuple val(meta), path("${meta.sample}_${meta.replicate}.bw"), emit: bigwig
    path("${meta.sample}_${meta.replicate}.bigwigtrack.txt"), emit: bigwigtrack

    script:
    """
    
    module load samtools
    module load bedtools
    module load ucscbigtools
    
    samtools sort -n -O BAM ${bam} > nsorted.bam
    samtools fixmate nsorted.bam matefixed.bam
    
    bedtools bamtobed -bedpe -i matefixed.bam > fragments.bedpe
    
    awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' fragments.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > fragments.bed
    
    bedtools genomecov -bg -i fragments.bed -g ${projectDir}/references/chromsizes/${meta.genome} > ${meta.sample}_${meta.replicate}.bg
    
    bedGraphToBigWig ${meta.sample}_${meta.replicate}.bg ${projectDir}/references/chromsizes/${meta.genome} ${meta.sample}_${meta.replicate}.bw
    
    printf 'track type=bigWig name=\"${meta.sample}_${meta.replicate}_bw\" bigDataUrl=cloudpath/${meta.sample}_${meta.replicate}.bw description=${meta.sample}_${meta.replicate} visibility=full smoothingWindow=4 windowingFunction=mean maxHeightPixels=100:50:8\\n' > ${meta.sample}_${meta.replicate}.bigwigtrack.txt
    
    """
}
