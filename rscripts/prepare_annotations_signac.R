# https://stuartlab.org/signac/articles/pbmc_vignette
# https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps

library(Signac)
library(AnnotationHub)

ah <- AnnotationHub()

# References 2024-A (March 13, 2024)
# Human reference, GRCh38 (GENCODE v44/Ensembl 110 annotations)
query(ah, "EnsDb.Hsapiens.v110")
ensdb_v110 <- ah[["AH113665"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v110)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
saveRDS(annotations, file = "hg38.rds")

# Mouse reference, GRCm39 (GENCODE v33/Ensembl 110 annotations)
query(ah, "EnsDb.Mmusculus.v110")
ensdb_v110 <- ah[["AH113713"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v110)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm39"
saveRDS(annotations, file = "mm39.rds")
