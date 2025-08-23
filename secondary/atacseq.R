library(magrittr)
library(tidyverse)
library(Signac)

setwd("~/GEDI/demo/atacseq")

WT_1 = read.table("WT_PP_1.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))

WT_1 %>% head()
WT_1 %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% head()
WT_1 %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% range()
WT_1 %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% mean()
WT_1 %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% median()

WT_2 = read.table("WT_PP_2.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))
TKO_1 = read.table("TKO_PP_1.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))
TKO_2 = read.table("TKO_PP_2.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))

####### Visualize #######
hg38_annotations = readRDS("hg38.rds")

genome(atac_small) = "hg38"
Annotation(atac_small) = hg38_annotations

WT_1_peaks = WT_1 %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()
WT_2_peaks = WT_2 %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()
TKO_1_peaks = TKO_1 %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()
TKO_2_peaks = TKO_2 %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()

my_window = "chr13-27913668-27932633"
gene_plot <- AnnotationPlot(atac_small, region = my_window)
WT_1_plot <- PeakPlot(atac_small, region = my_window, peaks = WT_1_peaks, color = "#00563E")
WT_2_plot <- PeakPlot(atac_small, region = my_window, peaks = WT_2_peaks, color = "#00563E")
TKO_1_plot <- PeakPlot(atac_small, region = my_window, peaks = TKO_1_peaks, color = "#860C1E")
TKO_2_plot <- PeakPlot(atac_small, region = my_window, peaks = TKO_2_peaks, color = "#860C1E")

CombineTracks(plotlist = list(gene_plot, WT_1_plot, WT_2_plot, TKO_1_plot, TKO_2_plot), heights = c(1, 1, 1, 1, 1))

WT_1_signal <- BigwigTrack(region = my_window, bigwig = "WT_PP_1.bw", y_label = "WT_1", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#00563E"))
WT_2_signal <- BigwigTrack(region = my_window, bigwig = "WT_PP_2.bw", y_label = "WT_2", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#00563E"))
TKO_1_signal <- BigwigTrack(region = my_window, bigwig = "TKO_PP_1.bw", y_label = "TKO_1", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#860C1E"))
TKO_2_signal <- BigwigTrack(region = my_window, bigwig = "TKO_PP_2.bw", y_label = "TKO_2", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#860C1E"))

CombineTracks(plotlist = list(gene_plot, WT_1_signal, WT_1_plot, WT_2_signal, WT_2_plot, TKO_1_signal, TKO_1_plot, TKO_2_signal, TKO_2_plot), 
              heights = c(1, 1, 1, 1, 1, 1, 1, 1, 1))

### Create a consensus (MERGE) - Continued in atacseq.sh ...

###
###
###

### ... continuing from atacseq.sh

merged = read.table("merged.bed") %>% 
  set_colnames(c("chr", "start", "end"))

merged %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% range()
merged %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% mean()
merged %>% mutate(peaksize = end-start) %>% pull(peaksize) %>% median()

####### Differential analysis (differential accessibility) #######

library(DESeq2)

counts = read.table("counts.tsv")

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

library(EnhancedVolcano)

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 1) %>% dim()
res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < -1) %>% dim()

# Save the peaks that are downregulated to make some heatmaps

ATAC_down_in_PP = res %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < -1) %>% 
  rownames_to_column("peak") %>% 
  separate(col = "peak", into = c("chr", "start", "end"), sep = "-") %>%
  arrange(chr, as.numeric(start)) %>%
  dplyr::select(chr, start, end)

write.table(ATAC_down_in_PP, "ATAC_down_in_PP.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

############# GREAT ###############

library(rGREAT)
set.seed(123)

dars = res %>% data.frame() %>% filter(padj < 0.05) %>% rownames() %>% StringToGRanges()
GREATres = great(dars, "GO:BP", "txdb:hg38")

plotRegionGeneAssociations(GREATres)

getEnrichmentTable(GREATres) %>% head()

############ Motif enrichments ??? ####











