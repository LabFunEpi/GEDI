#!/usr/bin/env nextflow

process BED_TO_BIGBED {

    publishDir "results/bigbeds", mode: 'copy', pattern: "*.bb"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bb"), emit: bigbed
    path("*.bigbedtrack.txt"), emit: bigbedtrack
    
    script:
    """

    if [ "${meta.replicate}" = "null" ]; then
        name=${meta.sample}
    else
        name=${meta.sample}_${meta.replicate}
    fi
    
    module load ucscbigtools
    
    awk -v OFS='\\t' '{print \$1,\$2,\$3}' ${bed} > temp.bed
    
    bedToBigBed temp.bed ${params.refDir}/chromsizes/${meta.genome} \${name}.bb
    
    printf 'track type=bigBed name=\"%s_bb\" bigDataUrl=cloudpath/%s.bb description=%s visibility=dense\\n' \\
        "\${name}" "\${name}" "\${name}" > \${name}.bigbedtrack.txt
    
    """
}
