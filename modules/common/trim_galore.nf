#!/usr/bin/env nextflow

process TRIM_GALORE {

    publishDir "results/qc/trimming_reports", mode: 'copy', pattern: "*_trimming_report.txt"
    publishDir "results/qc/fastqc_reports", mode: 'copy', pattern: "*_fastqc.{zip,html}"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sample}_${meta.replicate}_val_1.fq.gz"), path("${meta.sample}_${meta.replicate}_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
    path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    """
   
    module load fastqc
    module load cutadapt
   
    /research/bsi/tools/biotools/trim_galore/0.6.5/trim_galore \\
        --fastqc \\
        --paired ${fastq_1} ${fastq_2} \\
        --basename ${meta.sample}_${meta.replicate} \\
        --nextseq 20
   
    """
}
