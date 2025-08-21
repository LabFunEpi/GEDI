#!/usr/bin/env nextflow

process CREATE_SEURAT_MULT {

    publishDir "results", mode: 'copy', pattern: "${sample}_seurat.rds"

    input:
    tuple val(sample), val(genome), path(fbc_matrix), path(fragments), path(fragindex), path(metrics)

    output:
    path "${sample}_seurat.rds"

    script:
    """

    module load r

    Rscript ${projectDir}/rscripts/create_seurat_mult.R \\
        ${params.refDir} \\
        ${sample} \\
        ${genome} \\
        ${fbc_matrix} \\
        ${fragments} \\
        ${metrics} \\
        ${sample}_seurat.rds
    
    """
}
