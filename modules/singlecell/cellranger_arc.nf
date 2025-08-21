#!/usr/bin/env nextflow

process CELLRANGER_ARC {

    publishDir "results", mode: 'copy', pattern: "*/outs/*"
    publishDir "results", mode: 'copy', pattern: "${sample}.csv"

    input:
    tuple val(index), val(sample), val(genome)

    output:
    path "*/outs/*"
    path "${sample}.csv"
    tuple val(sample), val(genome), path("*/outs/filtered_feature_bc_matrix.h5"), path("*/outs/atac_fragments.tsv.gz"), path("*/outs/atac_fragments.tsv.gz.tbi"), path("*/outs/per_barcode_metrics.csv"), emit: toseurat

    script:
    """
    
    echo "fastqs,sample,library_type\n${projectDir}/data/gex,${sample},Gene Expression\n${projectDir}/data/atac,${sample},Chromatin Accessibility" \\
    > ${sample}.csv

    memory_with_unit="${task.memory}"
    memory_number=\$(echo "\$memory_with_unit" | sed 's/[^0-9.]//g')

    module load cellranger-arc
    
    cellranger-arc count \\
        --id=${sample} \\
        --reference=${params.refDir}/cellranger_arc/${genome} \\
        --libraries=${sample}.csv \\
        --localcores=$task.cpus \\
        --localmem=\$memory_number
    
    """
}
