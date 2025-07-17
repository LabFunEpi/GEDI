library("DESeq2")
setwd("~/GEDI")

### Read the data and define the sample conditions
directory <- "~/GEDI/results"

sampleFiles <- grep("htseq",list.files(directory),value=TRUE)

sampleNames <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6")
sampleConditions <- c("Control", "Control", "Control", "Treated", "Treated", "Treated")

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)

sampleTable$condition <- factor(sampleTable$condition, levels = c("Control", "Treated"))

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

### Pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

### Differential gene expression
dds <- DESeq(dds)
res <- results(dds)

### Shrink log fold change
resLFC <- lfcShrink(dds, coef="condition_Treated_vs_Control", type="apeglm")
plotMA(res, ylim=c(-6,6))
plotMA(resLFC, ylim=c(-6,6))

### Variance stabilization
vsd <- vst(dds, blind=FALSE)

library("vsn")
meanSdPlot(assay(normTransform(dds)))
meanSdPlot(assay(vsd))

### PCA
plotPCA(vsd, intgroup=c("condition"))

### Plot volcano plot
library(EnhancedVolcano)

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')




