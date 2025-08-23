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



