library(DESeq2)
library(tidyverse)

setwd("~/GEDI/demo/chipseq")

####### Differential analysis DE #######

counts = read.table("counts_DE.tsv")

counts_mtx = counts %>% dplyr::select(V4:V7) %>% as.matrix()
counts_mtx[is.na(counts_mtx)] = 0

peakNames = counts %>% dplyr::select(V1:V3) %>% unite(everything(), col = "peak", sep = "-") %>% pull(peak)

sampleNames = c("WT_DE_1", "WT_DE_2", "TKO_DE_1", "TKO_DE_2")
sampleConditions = c("WT", "WT", "TKO", "TKO")

sampleTable = data.frame(sampleName = sampleNames,
                         condition = sampleConditions)

sampleTable$condition = factor(sampleTable$condition, levels = c("WT", "TKO"))

colnames(counts_mtx) = sampleNames
rownames(counts_mtx) = peakNames

dds = DESeqDataSetFromMatrix(countData = round(counts_mtx),
                             colData = sampleTable,
                             design= ~ condition)

dds = DESeq(dds)
res = results(dds)

### Plot volcano plot
library(EnhancedVolcano)

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 1) %>% nrow()
res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < -1) %>% nrow()


###########################################################
###########################################################


####### Differential analysis PP #######

counts = read.table("counts_PP.tsv")

counts_mtx = counts %>% dplyr::select(V4:V7) %>% as.matrix()
counts_mtx[is.na(counts_mtx)] = 0

peakNames = counts %>% dplyr::select(V1:V3) %>% unite(everything(), col = "peak", sep = "-") %>% pull(peak)

sampleNames = c("WT_PP_1", "WT_PP_2", "TKO_PP_1", "TKO_PP_2")
sampleConditions = c("WT", "WT", "TKO", "TKO")

sampleTable = data.frame(sampleName = sampleNames,
                         condition = sampleConditions)

sampleTable$condition = factor(sampleTable$condition, levels = c("WT", "TKO"))

colnames(counts_mtx) = sampleNames
rownames(counts_mtx) = peakNames

dds = DESeqDataSetFromMatrix(countData = round(counts_mtx),
                             colData = sampleTable,
                             design= ~ condition)

dds = DESeq(dds)
res = results(dds)

### Plot volcano plot

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 1) %>% nrow()
res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < -1) %>% nrow()


###########################################################
###########################################################

# Save the peaks that are downregulated in FOXK2 binding to make some heatmaps

FOXA2_down_in_PP = res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < -1) %>% 
  rownames_to_column("peak") %>% 
  separate(col = "peak", into = c("chr", "start", "end"), sep = "-") %>%
  arrange(chr, as.numeric(start)) %>%
  dplyr::select(chr, start, end)

write.table(FOXA2_down_in_PP, "FOXA2_down_in_PP.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

### Continue in chipseq.sh ...

###
###
###

### Alternative to DESeq2, you could use DiffBind to do differential binding analysis
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html

### Coverage plot (Putting it all together)

library(Signac)

hg38_annotations = readRDS("~/GEDI/demo/atacseq/hg38.rds")

genome(atac_small) = "hg38"
Annotation(atac_small) = hg38_annotations

my_window = "chr7-127610289-127625877"
gene_plot <- AnnotationPlot(atac_small, region = my_window)

H3K4me1_WT_file = "H3K4me1_WT_PP_FC.bw"
H3K27ac_WT_file = "H3K27ac_WT_PP_FC.bw"
FOXA2_WT_file = "FOXA2_ChIP_WT_PP_FC.bw"
ATAC_WT_file = "~/GEDI/demo/atacseq/WT_PP_1.bw"
ATAC_TKO_file = "~/GEDI/demo/atacseq/TKO_PP_1.bw"
WGBS_WT_file = "~/GEDI/demo/dnamethyl/GSM4387092_WGBS_WT_PP.bw"
WGBS_TKO_file = "~/GEDI/demo/dnamethyl/GSM4387093_WGBS_TKO_PP.bw"

H3K4me1_WT <- BigwigTrack(region = my_window, bigwig = H3K4me1_WT_file, y_label = "H3K4me1", smooth = 250, downsample.rate = 1, ymax = 3) + 
  scale_fill_manual(values = c("#833a8c"))
H3K27ac_WT <- BigwigTrack(region = my_window, bigwig = H3K27ac_WT_file, y_label = "H3K27ac", smooth = 250, downsample.rate = 1, ymax = 20) + 
  scale_fill_manual(values = c("#5e5d35"))
FOXA2_WT <- BigwigTrack(region = my_window, bigwig = FOXA2_WT_file, y_label = "FOXA2", smooth = 250, downsample.rate = 1, ymax = 10) + 
  scale_fill_manual(values = c("#5f1319"))
ATAC <- BigwigTrack(region = my_window, bigwig = list(WT = ATAC_WT_file, TKO = ATAC_TKO_file), smooth = 250, downsample.rate = 1, bigwig.scale = "separate", ymax = 11) + 
  scale_fill_manual(values = c("#254a9b", "#254a9b"))
ATAC_WT <- BigwigTrack(region = my_window, bigwig = ATAC_WT_file, y_label = "ATAC_WT", smooth = 250, downsample.rate = 1, ymax = 6) + 
  scale_fill_manual(values = c("#254a9b"))
ATAC_TKO <- BigwigTrack(region = my_window, bigwig = ATAC_TKO_file, y_label = "ATAC_TKO", smooth = 250, downsample.rate = 1, ymax = 6) + 
  scale_fill_manual(values = c("#254a9b"))
WGBS_WT <- BigwigTrack(region = my_window, bigwig = WGBS_WT_file, y_label = "WGBS_WT", smooth = 250, downsample.rate = 1, ymax = 1) + 
  scale_fill_manual(values = c("black"))
WGBS_TKO <- BigwigTrack(region = my_window, bigwig = WGBS_TKO_file, y_label = "WGBS_TKO", smooth = 250, downsample.rate = 1, ymax = 1) + 
  scale_fill_manual(values = c("black"))

CombineTracks(plotlist = list(
  gene_plot, 
  H3K4me1_WT, 
  H3K27ac_WT, 
  FOXA2_WT, 
  ATAC, 
  WGBS_WT, 
  WGBS_TKO), 
  heights = c(1, 1, 1, 1, 2, 1, 1))






