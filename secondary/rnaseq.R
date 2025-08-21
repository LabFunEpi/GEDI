# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(DESeq2)
library(vsn)

setwd("~/GEDI/rnaseq")

### Read the data and define the sample conditions
directory = "~/GEDI/rnaseq"

sampleFiles = c("WT_ES_1_htseq.tab", "WT_ES_2_htseq.tab", 
                "WT_DE_1_htseq.tab", "WT_DE_2_htseq.tab", 
                "WT_GT_1_htseq.tab", "WT_GT_2_htseq.tab",
                "WT_PP_1_htseq.tab", "WT_PP_2_htseq.tab",
                "TKO2_ES_1_htseq.tab", "TKO2_ES_2_htseq.tab", 
                "TKO2_DE_1_htseq.tab", "TKO2_DE_2_htseq.tab", 
                "TKO2_GT_1_htseq.tab", "TKO2_GT_2_htseq.tab",
                "TKO2_PP_1_htseq.tab", "TKO2_PP_2_htseq.tab")

sampleNames = c("WT_ES_1", "WT_ES_2", 
                "WT_DE_1", "WT_DE_2", 
                "WT_GT_1", "WT_GT_2",
                "WT_PP_1", "WT_PP_2",
                "TKO2_ES_1", "TKO2_ES_2", 
                "TKO2_DE_1", "TKO2_DE_2", 
                "TKO2_GT_1", "TKO2_GT_2",
                "TKO2_PP_1", "TKO2_PP_2")

sampleConditions = c(rep("WT", 8), rep("TKO", 8))

sampleStages = rep(c("ES", "ES", "DE", "DE", "GT", "GT", "PP", "PP"), 2)

sampleTable = data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions,
                          stage = sampleStages)

sampleTable$condition = factor(sampleTable$condition, levels = c("WT", "TKO"))

dds = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

### Pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 1) >= smallestGroupSize
dds <- dds[keep,]

### Transformation
dds_t <- rlog(dds, blind=FALSE)

meanSdPlot(assay(normTransform(dds)))
meanSdPlot(assay(dds_t))

### PCA
plotPCA(dds_t, intgroup=c("condition", "stage"))

pcaData = plotPCA(dds_t, intgroup=c("condition", "stage"), returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

library(tidyverse)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=stage)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

### Differential gene expression

dds_PP = dds[, dds$stage %in% c("PP")]

dds_PP = DESeq(dds_PP)
resultsNames(dds_PP)

res = results(dds_PP)
resLFC = lfcShrink(dds_PP, coef = "condition_TKO_vs_WT", type="apeglm")

plotMA(res)
plotMA(resLFC)

### Volcano plot

library(EnhancedVolcano)

EnhancedVolcano(resLFC, lab = rownames(resLFC), x = 'log2FoldChange', y = 'pvalue')

resLFC %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange > 0) %>% nrow()
resLFC %>% data.frame() %>% filter(padj < 0.05 & log2FoldChange < 0) %>% nrow()

### Dealing with gene names

allgenes = rownames(dds)

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                  filters = "external_gene_name",
                  values = allgenes,
                  mart = mart)

listAttributes(mart)

### Pathway analysis (Pathway over-representation analysis)

# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# https://www.sciencedirect.com/science/article/pii/S2666675821000667

library(clusterProfiler)

deg = resLFC %>% data.frame() %>%
  filter(padj < 0.05) %>%
  rownames_to_column("external_gene_name") %>%
  left_join(gene_map) %>% 
  drop_na() %>%
  pull(entrezgene_id)

# search_kegg_organism('human', by='common_name')
gse <- enrichKEGG(deg, organism="hsa") %>% data.frame()

gseData = gse %>% head(n = 15) %>% dplyr::select(Description, p.adjust) %>% 
  mutate(negLogPadj = -log10(p.adjust))

ggplot(gseData, aes(x = negLogPadj, y = Description)) +
  geom_col()

ggplot(gseData, aes(x = negLogPadj, y = reorder(Description, negLogPadj))) +
  geom_col()

### Alternate - Pathway analysis (Gene set enrichment analysis - GSEA)

res_df = resLFC %>% data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(-log2FoldChange) %>%
  rownames_to_column("external_gene_name") %>%
  left_join(gene_map) %>% 
  drop_na()

sortedFoldchangeList = res_df$log2FoldChange
names(sortedFoldchangeList) = res_df$entrezgene_id

gse <- gseKEGG(geneList=sortedFoldchangeList, organism="hsa", by = "DOSE", pvalueCutoff = 1, nPerm = 100)

gse %>% data.frame() %>% head()

gseaplot(gse, geneSetID = 1, title = gse$Description[1])

### Alternate - Over-representation or GSEA on any geneset/collection of genesets (Generic)
# https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html

library(msigdbr)

msigdbr_species()
msigdbr_collections()

h_gene_sets = msigdbr(species = "Homo sapiens", collection = "H")
term2gene = h_gene_sets %>% dplyr::select(gs_name, gene_symbol)

# Over-representation analysis
deg = resLFC %>% data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

gse = enricher(deg, TERM2GENE = term2gene) %>% data.frame()

# GSEA
sortedFoldchangeList = res_df$log2FoldChange
names(sortedFoldchangeList) = res_df$external_gene_name

gse <- GSEA(sortedFoldchangeList, TERM2GENE = term2gene, by = "DOSE", pvalueCutoff = 1, nPerm = 100)

gse %>% data.frame() %>% head()

gseaplot(gse, geneSetID = 4, title = gse$Description[4])




