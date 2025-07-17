#!/usr/bin/env nextflow

process MACS2_CALLPEAKS {

    publishDir "results/beds", mode: 'copy', pattern: "${meta.sample}.macs2.bed"
    publishDir "results/qc/macs2_reports", mode: 'copy', pattern: "${meta.sample}.macs2.log"

    input:
    tuple val(meta), path(samplebam), path(controlbam)

    output:
    tuple val(meta), path("${meta.sample}.macs2.bed"), emit: macs2bed
    path "${meta.sample}.macs2.log", emit: macs2log

    script:
    """

    if [ "${meta.genome}" = "hg19" ] || [ "${meta.genome}" = "hg38" ]; then
        GSIZE="hs"
    else
        GSIZE="mm"
    fi
    
    /research/bsi/tools/biotools/macs2/2.2.9.1/bin/macs2 callpeak -t ${samplebam} -c ${controlbam} -f BAM -g \${GSIZE} -n ${meta.sample} --keep-dup all --nomodel 2> ${meta.sample}.macs2.log
    
    awk '{print \$1,\$2,\$3,\$4}' OFS='\\t' ${meta.sample}_peaks.narrowPeak > ${meta.sample}.macs2.bed
    
    """
}
