#!/usr/bin/env nextflow

include { CELLRANGER_ARC } from './modules/singlecell/cellranger_arc.nf'
include { CREATE_SEURAT_MULT } from './modules/singlecell/create_seurat_mult.nf'

// Primary input
params.input_csv = "data/samples.csv"

workflow {

    // Read and parse input
    def index = 0
    samplesheet_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row ->
            [index:++index, sample:row.sample, genome:row.genome]
        }
    
    // Main pipeline
    CELLRANGER_ARC(samplesheet_ch)
    CREATE_SEURAT_MULT(CELLRANGER_ARC.out.toseurat)

}
