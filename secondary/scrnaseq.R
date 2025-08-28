# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217837

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(BuenColors)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

#################### MY FIRST SEURAT OBJECT ##########################

setwd("~/GEDI/demo/scrnaseq")

Sys.setenv(VROOM_CONNECTION_SIZE=500000)
scRNA_df <- readr::read_csv("GSE217837_scRNA_count_matrix.csv.gz", col_names = TRUE)
scRNA_df <- column_to_rownames(scRNA_df, "...1")

# Or load data from alternate formats
# scRNA_df <- Read10X_h5("filtered_feature_bc_matrix.h5")
# scRNA_df <- Read10X("filtered_feature_bc_matrix")


scRNA_df <- Matrix::Matrix(as.matrix(scRNA_df),sparse = T)

scRNA_sobj <- CreateSeuratObject(counts = scRNA_df)

VlnPlot(scRNA_sobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

scRNA_sobj <- NormalizeData(scRNA_sobj)
scRNA_sobj <- FindVariableFeatures(scRNA_sobj)
scRNA_sobj <- ScaleData(scRNA_sobj)
scRNA_sobj <- RunPCA(scRNA_sobj)
scRNA_sobj <- RunUMAP(scRNA_sobj, dims = 1:50)
scRNA_sobj <- FindNeighbors(scRNA_sobj)
scRNA_sobj <- FindClusters(scRNA_sobj)

DimPlot(scRNA_sobj, reduction = "umap")

# AvailableData()$Dataset
scRNA_sobj <- RunAzimuth(scRNA_sobj, reference = "pancreasref")

DimPlot(scRNA_sobj, group.by = "predicted.annotation.l1") + 
  scale_color_manual(values = jdb_palette("corona"))

Idents(scRNA_sobj) <- "predicted.annotation.l1"
VlnPlot(scRNA_sobj, features = c("FOXA2", "PAX4"))

#################### MY SECOND SEURAT OBJECT ##########################

snRNA_df <- readr::read_csv("GSE217837_snRNA_count_matrix.csv.gz", col_names = TRUE)
snRNA_df <- column_to_rownames(snRNA_df, "...1")
snRNA_df <- Matrix::Matrix(as.matrix(snRNA_df),sparse = T)
snRNA_sobj <- CreateSeuratObject(counts = snRNA_df)

snRNA_sobj <- snRNA_sobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50) %>%
  FindNeighbors() %>% 
  FindClusters() %>% 
  RunAzimuth(reference = "pancreasref")

DimPlot(snRNA_sobj, group.by = "predicted.annotation.l1") + 
  scale_color_manual(values = jdb_palette("corona"))

Idents(snRNA_sobj) <- "predicted.annotation.l1"
VlnPlot(snRNA_sobj, features = c("FOXA2", "PAX4"))

#################### PUTTING THEM TOGETHER - BATCH CORRECTION ##########################

scRNA_sobj$group <- "sc"
snRNA_sobj$group <- "sn"
combined_sobj <- merge(
  x = scRNA_sobj, 
  y = snRNA_sobj, 
  add.cell.ids = c("sc", "sn")
)

combined_sobj <- combined_sobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = "umap.uncorrected") %>% 
  FindNeighbors() %>% 
  FindClusters()

DimPlot(combined_sobj, reduction = "umap.uncorrected", group.by = c("group", "predicted.annotation.l1")) + 
  scale_color_manual(values = jdb_palette("corona"))

combined_sobj <- IntegrateLayers(
  object = combined_sobj, 
  method = CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "cca"
)
combined_sobj <- RunUMAP(
  object = combined_sobj, 
  reduction = "cca", 
  dims = 1:50, 
  reduction.name = "umap.corrected"
)
combined_sobj <- JoinLayers(combined_sobj)

DimPlot(combined_sobj, reduction = "umap.corrected", group.by = c("group", "predicted.annotation.l1")) + 
  scale_color_manual(values = jdb_palette("corona"))

#################### DIFFERENTIAL GENE EXPRESSION ##########################

Idents(combined_sobj) <- paste(combined_sobj$group, combined_sobj$predicted.annotation.l1, sep = ".")
degs <- FindMarkers(
  object = combined_sobj, 
  ident.1 = "sn.beta", 
  ident.2 = "sc.beta",
  min.pct = 0.4
)
head(degs)

ggplot(data = degs) +
  geom_point(mapping = aes(x = avg_log2FC, y = -log10(p_val_adj)))

######################### RESOURCES LINKS ###############################

# https://www.10xgenomics.com/blog/single-cell-rna-seq-an-introductory-overview-and-tools-for-getting-started
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
# https://www.10xgenomics.com/products/loupe-browser
# https://scanpy.readthedocs.io/en/stable/
# https://satijalab.org/seurat/
# https://satijalab.org/seurat/articles/get_started_v5
# https://github.com/swolock/scrublet
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
# https://azimuth.hubmapconsortium.org/
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://github.com/dviraran/SingleR
# https://github.com/PaulingLiu/scibet
# http://velocyto.org/
# https://github.com/sqjin/CellChat
# https://cole-trapnell-lab.github.io/cicero-release/docs_m3/
# https://github.com/GreenleafLab/chromVAR
# https://github.com/buenrostrolab/FigR
# https://singlecell.broadinstitute.org/single_cell
# https://www.ebi.ac.uk/gxa/sc/home


# Single cell best practices article: https://www.nature.com/articles/s41576-023-00586-w
# https://www.sc-best-practices.org/

# RStudio in the cloud
# https://posit.cloud/

# Learn bioinformatics through Linux cloud
# https://sandbox.bio/

# scRNA_df <- Read10X_h5("filtered_feature_bc_matrix.h5")
# scRNA_df <- Read10X("filtered_feature_bc_matrix")
