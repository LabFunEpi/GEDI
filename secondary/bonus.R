
###### GSEA alternate #################
# https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html

DEG = read.csv("DEG.csv")
res_df = DEG %>% left_join(gene_entrez_map, by = c("Gene" = "SYMBOL")) %>% drop_na() %>%
  mutate(FC = 2^log2FC) %>%
  arrange(-log2FC)

geneList = res_df$log2FC
names(geneList) = res_df$ENTREZID

gse <- enrichKEGG(names(geneList), organism="hsa")

# gse <- gseKEGG(geneList=geneList, organism="hsa", by = "DOSE", pvalueCutoff = 1, nPerm = 100)
# data.frame(gse) %>% arrange(pvalue) %>% head(n = 10)
# data.frame(gse) %>% dim()
# data.frame(gse) %>% filter(ID %in% c("hsa04911", "hsa04950"))

# KEGG = download_KEGG(species="hsa")
# "hsa04911", "hsa04950"

res_df = res %>% data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(-log2FoldChange) %>%
  rownames_to_column("SYMBOL") %>%
  left_join(gene_entrez_map) %>% 
  drop_na()

geneList = res_df$log2FoldChange
names(geneList) = res_df$ENTREZID

gse <- gseKEGG(geneList=geneList, organism="hsa", by = "DOSE", pvalueCutoff = 1, nPerm = 100)

####### Visualize using Gviz #######
library(Gviz)
library(GenomicRanges)
library(rtracklayer)

WT_bed = AnnotationTrack("WT_PP.macs2.bed", name = "WT_PP", 
                         genome = "hg38", chromosome = "chr13", 
                         start = 27913668, end = 27932633)
TKO_bed = AnnotationTrack("TKO_PP.macs2.bed", name = "TKO_PP", 
                          genome = "hg38", chromosome = "chr13", 
                          start = 27913668, end = 27932633)


grtrack = GeneRegionTrack("hg38.gtf", 
                          name = "Genes", genome = "hg38", 
                          chromosome = "chr13", start = 27913668, 
                          end = 27932633, transcriptAnnotation = "gene_name", 
                          stacking = 'squish', stackHeight = 0.3)

gtrack = GenomeAxisTrack(genome = "hg38", chromosome = "chr13", 
                         start = 27913668, end = 27932633)


plotTracks(list(grtrack, WT_bed, TKO_bed, gtrack), chromosome = "chr13",
           from = 27913668, to = 27932633)

WT_bw_1 = DataTrack("WT_PP_1.bw", type = "histogram", 
                    name = "WT_PP_1.bw", genome = "hg38", chromosome = "chr13", 
                    start = 27913668, end = 27932633)
WT_bw_2 = DataTrack("WT_PP_2.bw", type = "histogram", 
                    name = "WT_PP_2.bw", genome = "hg38", chromosome = "chr13", 
                    start = 27913668, end = 27932633)
TKO_bw_1 = DataTrack("TKO_PP_1.bw", type = "histogram", 
                     name = "TKO_PP_1.bw", genome = "hg38", chromosome = "chr13", 
                     start = 27913668, end = 27932633)
TKO_bw_2 = DataTrack("TKO_PP_2.bw", type = "histogram", 
                     name = "TKO_PP_2.bw", genome = "hg38", chromosome = "chr13", 
                     start = 27913668, end = 27932633)

plotTracks(list(grtrack, WT_bw_1, WT_bw_2, WT_bed, TKO_bw_1, TKO_bw_2, TKO_bed, 
                gtrack), chromosome = "chr13", from = 27913668, to = 27932633)

############# Visualize averaged bigwigs ###########

WT = read.table("WT_PP.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))

TKO = read.table("TKO_PP.macs2.bed") %>% 
  set_colnames(c("chr", "start", "end", "id"))

WT_peaks = WT %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()
TKO_peaks = TKO %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()
merged_peaks = merged %>% unite("peaks", chr:end, sep = "-") %>% pull(peaks) %>% StringToGRanges()

WT_plot <- PeakPlot(atac_small, region = my_window, peaks = WT_peaks, color = "#00563E")
TKO_plot <- PeakPlot(atac_small, region = my_window, peaks = TKO_peaks, color = "#00563E")
merged_plot <- PeakPlot(atac_small, region = my_window, peaks = merged_peaks, color = "black")

WT_signal <- BigwigTrack(region = my_window, bigwig = "WT_PP.bw", y_label = "WT", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#00563E"))
TKO_signal <- BigwigTrack(region = my_window, bigwig = "TKO_PP.bw", y_label = "TKO", downsample.rate = 1, bigwig.scale = "separate") + 
  scale_fill_manual(values = c("#860C1E"))

CombineTracks(plotlist = list(gene_plot, WT_signal, WT_plot, TKO_signal, TKO_plot, merged_plot), heights = c(1, 1, 1, 1, 1, 1))

