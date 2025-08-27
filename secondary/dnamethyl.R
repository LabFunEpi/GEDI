
# https://groups.google.com/g/methylkit_discussion/c/7cX28UgQ0js?pli=1
# remotes::install_version("data.table", version = "1.16")
library(data.table, lib.loc="~/R/x86_64-pc-linux-gnu-library/4.4")

packageVersion("data.table")
library(methylKit)

library(magrittr)
library(tidyverse)
library(ggrastr)

setwd("~/GEDI/demo/dnamethyl")

file.list = list("SRR11248464_chr7_CpG.methylKit", 
                 "SRR11248465_chr7_CpG.methylKit")

myobj=methRead(file.list,
               sample.id=list("WT_PP","TKO_PP"),
               assembly="hg38",
               treatment=c(0,1),
               mincov = 10)

getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

meth=methylKit::unite(myobj, destrand=FALSE)

getCorrelation(meth,plot=TRUE)

myDiff=calculateDiffMeth(meth) 

myDiff = myDiff %>% data.frame() %>%
  mutate(negLogQvalue = -log10(qvalue)) %>%
  unite(chr:start, col = "CpGsite", sep = "-")

### Plot volcano plot

ggplot(myDiff, aes(x = meth.diff, y = negLogQvalue)) +
  geom_point_rast(data = myDiff %>% filter(abs(meth.diff) < 50), size=1, stroke=0.1, shape=16) +
  geom_point_rast(data = myDiff %>% filter(meth.diff > 50), size=1, stroke=0.1, shape=16, color = "red") +
  geom_point_rast(data = myDiff %>% filter(meth.diff < -50), size=1, stroke=0.1, shape=16, color = "blue") +
  theme_classic()

  
myDiff %>% filter(pvalue < 0.05 & meth.diff > 50) %>% dim()
myDiff %>% filter(pvalue < 0.05 & meth.diff < -50) %>% dim()


### ReMapEnrich analysis 

library(devtools)
install_github("remap-cisreg/ReMapEnrich")
library(ReMapEnrich)

setwd("~/GEDI/demo/dnamethyl")

query <- bedToGranges("~/GEDI/demo/dnamethyl/Ex2_GEDI_HyperDMR.bed")
remapCatalog <- bedToGranges("~/GEDI/demo/dnamethyl/Ex2_remap2022_pancreatic-progenitor_nr_macs2_hg38_v1_0.bed")

enrichment.df <- enrichment(query, remapCatalog)

enrichmentDotPlot(enrichment.df)

