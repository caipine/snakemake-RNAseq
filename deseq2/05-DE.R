load("results/preprocessing_2.RData")
load("results/annot.RData")
library(dplyr)
library(biomaRt)
library(DESeq2)
library(tidyverse)


resultsNames(ddsObj)
ddsObj$condition                             

ddsObj$condition <- relevel(ddsObj$condition, ref = "Sensitive")

resultsNames(ddsObj)
ddsObj$condition


ddsShrink_Condition_tmp <- lfcShrink(ddsObj, coef="condition_Resistant_vs_Sensitive")
#or load("results/ddsShrink_Condition_tmp.RData")
save(ddsShrink_Condition_tmp, file="results/ddsShrink_Condition_tmp.RData")
load("results/ddsShrink_Condition_tmp.RData")


shrinkLvV_Condition <- as.data.frame(ddsShrink_Condition_tmp) %>%
    rownames_to_column("GeneID") %>% 
    left_join(annotL, "GeneID") 
dim(shrinkLvV_Condition)
shrinkLvV_Condition <- shrinkLvV_Condition[!is.na(shrinkLvV_Condition$Symbol),]
dim(shrinkLvV_Condition)

library("ggrepel")
#shrinkLvV_Condition$Significant <- ifelse(shrinkLvV_Condition$padj < 0.05, "FDR < 0.1", "Not Sig")
#shrinkLvV_Condition$Significant <- ifelse((shrinkLvV_Condition$padj < 0.05 & shrinkLvV_Condition$padj  > 0.01) & abs(shrinkLvV_Condition$log2FoldChange) >1 , "0.01 < FDR < 0.05", "Not_Sig")

shrinkLvV_Condition <- shrinkLvV_Condition[!is.na(shrinkLvV_Condition$padj),]

shrinkLvV_Condition$Significant <- ifelse((shrinkLvV_Condition$padj < 0.05 & shrinkLvV_Condition$padj  > 0.01)  , "0.01 < FDR < 0.05", "Not_Sig")
shrinkLvV_Condition$Significant2 <- ifelse((shrinkLvV_Condition$padj <= 0.01 )  , "FDR < 0.01", "Not_Sig")

shrinkLvV_Condition$Significant[shrinkLvV_Condition$Significant2 == "FDR < 0.01"] <- shrinkLvV_Condition$Significant2[shrinkLvV_Condition$Significant2 == "FDR < 0.01"] 

head(shrinkLvV_Condition)

dim(shrinkLvV_Condition[shrinkLvV_Condition$padj <= 0.01,])
#[1] 567  12

dim(shrinkLvV_Condition[shrinkLvV_Condition$padj > 0.01 & shrinkLvV_Condition$padj <= 0.05,])
#[1] 1189  12


dim(shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.01 & abs(shrinkLvV_Condition$log2FoldChange) >1,])
#[1]497 12



dim(shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.01 & shrinkLvV_Condition$log2FoldChange < -1,])
#[1]344 12
shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.01 & shrinkLvV_Condition$log2FoldChange < -1,]$Symbol


unique(shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.01 & shrinkLvV_Condition$log2FoldChange < -1,]$Symbol)

#301
up <- unique(shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.0001 & shrinkLvV_Condition$log2FoldChange > 1,]$Symbol)
write.table(up, file = "up301.txt", row.names =F, quote =F, col.names=F)




tiff("IBN00~condition_DE.tif",width=1000,height=1000)
ggplot(shrinkLvV_Condition, aes(x = log2FoldChange, y= -log10(padj) )) + 
    geom_point(aes(colour= Significant), shape=20, size=3) +
    scale_color_manual(values = c("green","red","grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_hline(yintercept=2, colour="blue", linetype="dashed") +
    geom_vline(xintercept = -1:1, colour="red", linetype="dashed")
dev.off()





tiff("IBN00~condition_DE_label0.05_2.tif",width=1000,height=1000)
ggplot(shrinkLvV_Condition, aes(x = log2FoldChange, y= -log10(padj) )) + 
    geom_point(aes(colour= Significant), shape=20, size=3) +
    scale_color_manual(values = c("green","red","grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_hline(yintercept=2, colour="blue", linetype="dashed") +
    geom_vline(xintercept = -1:1, colour="red", linetype="dashed")+
   geom_text_repel(
    data = subset(shrinkLvV_Condition, (padj < 0.0001 & abs(log2FoldChange) >1 )),
    aes(label = Symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
)
dev.off()




pdf("IBN00-Condition_DE_label2.pdf",width=20,height=20)
ggplot(shrinkLvV_Condition, aes(x = log2FoldChange, y= -log10(padj) )) + 
    geom_point(aes(colour= Significant), shape=20, size=3) +
    scale_color_manual(values = c("red", "grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_hline(yintercept=2, colour="blue", linetype="dashed") +
    geom_vline(xintercept = -1:1, colour="red", linetype="dashed")+
   geom_text_repel(
    data = subset(shrinkLvV_Condition, (padj < 0.01 )),
    aes(label = Symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
)
dev.off()


