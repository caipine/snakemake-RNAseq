load("results/preprocessing_2.RData")
load("results/annot.RData")
library(dplyr)
#library(biomaRt)
library(DESeq2)
library(tidyverse)
library(gplots)
library(RColorBrewer)

resCond_vV <- results(ddsObj, name="condition_Resistant_vs_Sensitive")
#load("results/resCond_vV.data")

resCond_vV  <- as.data.frame(resCond_vV) %>% 
    rownames_to_column("GeneID") %>% 
    left_join( annotL, "GeneID") 

hmDat <- logNormalizedCounts[resCond_vV$GeneID,]
hmDat <- data.frame(round(hmDat, 2))
hmDat$NAME <- resCond_vV$Symbol
hmDat$Description <- resCond_vV$description
hmDat <- hmDat[,c(17:18,1:16)]
head(hmDat)
hmDat <- hmDat[!is.na(hmDat$NAME),]



write.table("#1.2", "IBN00.gct", quote =F, row.names = F, col.names=F)
l2 <- t(data.frame( c(nrow(hmDat), ncol(hmDat)-2)))
write.table(l2, "IBN00.gct",append = T,  quote =F, row.names = F, col.names=F, sep=" ")
write.table(hmDat, "IBN00.gct", append = T, quote =F, row.names = F, sep= "\t")


cls_l2 <- t(data.frame( c(ncol(hmDat)-2, 2, 1)))
write.table(cls_l2, "IBN00.cls", quote =F, row.names = F, col.names=F)
write.table("#Resistant Sensitive", "IBN00.cls",append = T,  quote =F, row.names = F, col.names=F, sep=" ")
write.table( t(data.frame(sampleinfo$condition)), "IBN00.cls", append = T, sep=" ", row.names=F,col.names=F, quote = F)

