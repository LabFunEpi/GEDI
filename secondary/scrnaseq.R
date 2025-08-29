# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217837

devtools::install_github("caleblareau/BuenColors")
devtools::install_github("chloelulu/scMayoMap")

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(EnhancedVolcano)
library(BuenColors)
library(scMayoMap)

options(future.globals.maxSize = 1e9)
set.seed(123)

setwd("~/GEDI/demo/scrnaseq")

# Read the data
pbmc.data = Read10X(data.dir = "~/GEDI/demo/scrnaseq/filtered_gene_bc_matrices/hg19")

# Or load data from alternate formats
# pbmc.data = Read10X_h5("filtered_feature_bc_matrix.h5")

# Create Seurat object
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Plot some QC
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Slicing the counts matrix

pbmc$RNA$counts[1:20, 1:20]

pbmc$RNA$counts[c("GAPDH", "ACTB"), 1:20]

pbmc$RNA$counts["GAPDH", ]

pbmc$RNA$counts[, 1]


# Why do we need to do Log Normalization?

hist(pbmc$RNA$counts["GAPDH",], breaks = 100)

pbmc = NormalizeData(pbmc)

hist(pbmc$RNA$data["GAPDH",], breaks = 100)


# Scaling the data: mean 0 and variance 1; why do we do this?

pbmc = ScaleData(pbmc, features = rownames(pbmc))

hist(pbmc$RNA$scale.data["GAPDH",], breaks = 100)


# Why do we want the most variable genes?

pbmc = FindVariableFeatures(pbmc)

VariableFeaturePlot(pbmc)


# PCA (First dimensionality reduction: Linear)

pbmc = RunPCA(pbmc)

DimPlot(pbmc, reduction = "pca") + NoLegend()

# temp = RunPCA(pbmc, features = rownames(pbmc))
# DimPlot(temp, reduction = "pca") + NoLegend()

ElbowPlot(pbmc)


# UMAP (Second dimensionality reduction: Non-linear)

pbmc = RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")


# Clustering

pbmc = FindNeighbors(pbmc, dims = 1:10)

pbmc = FindClusters(pbmc, resolution = 0.5)

DimPlot(pbmc, reduction = "umap")


# Find markers for the clusters

DefaultAssay(pbmc) = "RNA"
Idents(pbmc) = "seurat_clusters"

markers = FindAllMarkers(pbmc, method = 'MAST')

topmarkers = markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) %>% 
  pull(gene) %>% 
  unique()

DotPlot(pbmc, features = topmarkers) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(pbmc, features = "CD79A")


# Identify cell-types that match the markers

scMayoMap.obj = scMayoMap(data = markers, database=scMayoMapDatabase, tissue = 'blood')

scMayoMap.plot(scMayoMap.object = scMayoMap.obj)


# Assign cell-types to the clusters

pbmc$celltype = case_match(
  pbmc$seurat_clusters,
  "0" ~ "CD4 Naive T cell", 
  "1" ~ "CD14 Monocyte", 
  "2" ~ "CD4 Central Memory T cell", 
  "3" ~ "Naive B cell", 
  "4" ~ "CD8 Effector Memory T cell",
  "5" ~ "CD16 Monocyte",
  "6" ~ "CD56-dim natural killer cell", 
  "7" ~ "Plasmacytoid dendritic cell",
  "8" ~ "Megakaryocyte"
)

Idents(pbmc) = "celltype"

DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()


# Visualize some markers

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

VlnPlot(pbmc, features = c("MS4A1", "GNLY"))

topmarkers = markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) %>% 
  pull(gene) %>% 
  unique()

DotPlot(pbmc, features = topmarkers) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Differential gene expression

degs = FindMarkers(
  object = pbmc, 
  ident.1 = "CD14 Monocyte",
  ident.2 = "CD16 Monocyte",
  min.pct = 0.4
)

head(degs)

EnhancedVolcano(degs, lab = rownames(degs), x = 'avg_log2FC', y = 'p_val_adj')


#################### PANCREAS DATA - OBJECT 1 ##########################

setwd("~/GEDI/demo/scrnaseq")

Sys.setenv(VROOM_CONNECTION_SIZE=500000)
scRNA_df = readr::read_csv("GSE217837_scRNA_count_matrix.csv.gz", col_names = TRUE)
scRNA_df = column_to_rownames(scRNA_df, "...1")
scRNA_df = Matrix::Matrix(as.matrix(scRNA_df),sparse = T)

scRNA_sobj = CreateSeuratObject(counts = scRNA_df)

scRNA_sobj = scRNA_sobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>%
  FindNeighbors() %>% 
  FindClusters(resolution = 0.2)

markers = FindAllMarkers(scRNA_sobj, method = 'MAST')
scMayoMap.obj = scMayoMap(data = markers, database=scMayoMapDatabase, tissue = 'pancreas')
scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

DimPlot(scRNA_sobj) + scale_color_manual(values = jdb_palette("corona"))

#################### PANCREAS DATA - OBJECT 2 ##########################

snRNA_df = readr::read_csv("GSE217837_snRNA_count_matrix.csv.gz", col_names = TRUE)
snRNA_df = column_to_rownames(snRNA_df, "...1")
snRNA_df = Matrix::Matrix(as.matrix(snRNA_df),sparse = T)

snRNA_sobj = CreateSeuratObject(counts = snRNA_df)

snRNA_sobj = snRNA_sobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>%
  FindNeighbors() %>% 
  FindClusters(resolution = 0.2)

markers = FindAllMarkers(snRNA_sobj, method = 'MAST')
scMayoMap.obj = scMayoMap(data = markers, database=scMayoMapDatabase, tissue = 'pancreas')
scMayoMap.plot(scMayoMap.object = scMayoMap.obj)

DimPlot(snRNA_sobj) + scale_color_manual(values = jdb_palette("corona"))

#################### PUTTING THEM TOGETHER - BATCH CORRECTION ##########################
# https://satijalab.org/seurat/articles/seurat5_integration

scRNA_sobj$group = "sc"
snRNA_sobj$group = "sn"
combined_sobj = merge(
  x = scRNA_sobj, 
  y = snRNA_sobj, 
  add.cell.ids = c("sc", "sn")
)

combined_sobj = combined_sobj %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10, reduction.name = "umap.uncorrected")

DimPlot(combined_sobj, reduction = "umap.uncorrected", group.by = "group") + scale_color_manual(values = jdb_palette("corona"))

combined_sobj = IntegrateLayers(
  object = combined_sobj, 
  method = RPCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "rpca"
)
combined_sobj = RunUMAP(
  object = combined_sobj, 
  reduction = "rpca", 
  dims = 1:10, 
  reduction.name = "umap.corrected"
)
combined_sobj = JoinLayers(combined_sobj)

DimPlot(combined_sobj, reduction = "umap.corrected", group.by = "group") + 
  scale_color_manual(values = jdb_palette("corona"))


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

# scRNA_df = Read10X_h5("filtered_feature_bc_matrix.h5")
# scRNA_df = Read10X("filtered_feature_bc_matrix")
