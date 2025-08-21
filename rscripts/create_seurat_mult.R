library(Seurat)
library(Signac)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
refDir <- args[1]
sample <- args[2]
genome <- args[3]
fbc_matrix <- args[4]
fragments <- args[5]
metrics <- args[6]
outfile <- args[7]

inputdata.10x <- Read10X_h5(fbc_matrix)
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
metadata <- read.csv(
    file = metrics,
    header = TRUE,
    row.names = 1
)

sobj <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
if(genome == "hg38"){
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
}else{
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
}

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

annotations <- readRDS(paste0(refDir, "/annotations/", genome, ".rds"))

chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = genome,
    fragments = fragments,
    annotation = annotations
)
sobj[["ATAC"]] <- chrom_assay

sobj <- NucleosomeSignal(object = sobj, assay = "ATAC")
sobj <- TSSEnrichment(object = sobj, assay = "ATAC")

saveRDS(sobj, outfile)
