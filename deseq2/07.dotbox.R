load("results/preprocessing_2.RData")
load("results/annot.RData")
library(dplyr)
#library(biomaRt)
library(DESeq2)
library(tidyverse)
library(gplots)
library(RColorBrewer)


load("results/ddsShrink_Condition_tmp.RData")


shrinkLvV_Condition <- as.data.frame(ddsShrink_Condition_tmp) %>%
    rownames_to_column("GeneID") %>% 
    left_join(annotL, "GeneID") 
dim(shrinkLvV_Condition)
shrinkLvV_Condition <- shrinkLvV_Condition[!is.na(shrinkLvV_Condition$Symbol),]
dim(shrinkLvV_Condition)

up <- unique(shrinkLvV_Condition[shrinkLvV_Condition$padj <  0.0001 & shrinkLvV_Condition$log2FoldChange > 1,]$Symbol)
up <- up[!is.na(up)]



resCond_vV <- results(ddsObj, name="condition_Resistant_vs_Sensitive")
resCond_vV <-  data.frame(resCond_vV) %>% 
    rownames_to_column("GeneID") %>% 
    left_join( annotL, "GeneID") 




normalizedCounts <- counts(ddsObj, normalized=TRUE)
hmCount  <- as.data.frame(normalizedCounts) %>% 
    rownames_to_column("GeneID") %>% 
    left_join( annotL, "GeneID") 



###################
#box dot


namei <- "DGAT2"
t1 <- resCond_vV[resCond_vV$Symbol %in% namei,] 
t1_geneid<- t1[order(t1$pvalue),][1,1]
t2 <- data.frame(normalizedCounts[rownames(normalizedCounts) %in% t1_geneid,])
colnames(t2) <- "reads"
t2$sample <- rownames(t2)
t2$gene <- namei
t2$condition <- sampleinfo$condition
#t2
t3 <- t2



tiff("IBN00~condition_DGAT2.boxdot.tif",width=1000,height=1000)
ggplot( t3, aes(x = condition, y = reads ) ) + 
  geom_boxplot(    aes(color=condition, fill=condition), 
    alpha=.25, outlier.size=1, outlier.colour="gray" ) + 
  ylab("reads") + 
      		geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="black")+
  facet_wrap( ~ gene , scales = "free_y")
dev.off()



for (i in up) {
namei <-i
t1 <- resCond_vV[resCond_vV$Symbol %in% namei,] 
t1_geneid<- t1[order(t1$pvalue),][1,1]
t2 <- data.frame(normalizedCounts[rownames(normalizedCounts) %in% t1_geneid,])
colnames(t2) <- "reads"
t2$sample <- rownames(t2)
t2$gene <- namei
t2$condition <- sampleinfo$condition
t3 <- rbind(t3,t2)
}


tiff("IBN00~condition_gene.boxdot.tif",width=2000,height=1000)
ggplot( t3[18:nrow(t3),], aes(x = condition, y = reads ) ) + 
	geom_boxplot(    aes(color=condition, fill=condition), 
    		alpha=.25, outlier.size=1, outlier.colour="gray" ) + 
	ylab("reads") + 
      	geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="black")+
	theme(strip.text = element_text(size=25)) +
	facet_wrap( ~ gene , scales = "free_y", nrow = NULL, ncol = 9)
dev.off()


tiff("IBN00~condition_geneEv.tif",width=2000,height=1000)
ggplot( t3[18:nrow(t3),], aes(x = condition, y = reads ) ) + 
	geom_violin(trim = FALSE, aes(color=condition, fill=condition), alpha=.25, outlier.size=1, outlier.colour="gray")+
	ylab("reads") + 
	geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="black")+
	theme(strip.text = element_text(size=25)) +
	facet_wrap( ~ gene , scales = "free_y", ncol = 9)

dev.off()




