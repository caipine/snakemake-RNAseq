load("results/preprocessing_2.RData")
load("results/annot.RData")
library(dplyr)
#library(biomaRt)
library(DESeq2)
library(tidyverse)
library(gplots)
library(RColorBrewer)

resultsNames(ddsObj)
ddsObj$condition                             

ddsObj$condition <- relevel(ddsObj$condition, ref = "Sensitive")

resultsNames(ddsObj)
ddsObj$condition





resCond_vV <- results(ddsObj, name="condition_Resistant_vs_Sensitive")
#load("results/resCond_vV.data")

resCond_vV  <- as.data.frame(resCond_vV) %>% 
    rownames_to_column("GeneID") %>% 
    left_join( annotL, "GeneID") 




resCond_vV.nona<- resCond_vV[!is.na(resCond_vV$Symbol),]
resCond_vV.0.01 <- resCond_vV.nona[(resCond_vV.nona$padj < 0.01),]
resCond_vV.0.01 <- resCond_vV.0.01[!is.na(resCond_vV.0.01$GeneID),]
dim(resCond_vV.0.01) 

length(unique(resCond_vV.0.01$Symbol)) # 463

length(unique(resCond_vV.0.01[resCond_vV.0.01$log2FoldChange > 0,]$Symbol))  # up 314

length(unique(resCond_vV.0.01[resCond_vV.0.01$log2FoldChange < 0,]$Symbol))  # down 149



hmDat <- logNormalizedCounts[resCond_vV.0.01$GeneID,]

# Get some nicer colours
mypalette <- brewer.pal(11, "RdYlBu")
# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.sample <- c("purple","orange")[sampleinfo$condition]
# Plot the heatmap

# Get some nicer colours
mypalette <- brewer.pal(11, "RdYlBu")
# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.sample <- c("purple","orange")[sampleinfo$condition]
# Plot the heatmap



condition_color <- unlist( lapply(sampleinfo$condition, function(x){  
                  if (grepl("Resistant", x)) "purple"     
                  else if (grepl("Sensitive", x)) "orange"
                  }))

tissue_color <- unlist( lapply(sampleinfo$tissue, function(x){  
                  if (grepl("Spleen&LN", x)) "red"     
                  else if (grepl("PB", x)) "green"
                  else if (grepl("Apheresis", x)) "yellow"
                  else if (grepl("Spleen", x)) "blue"
                  }))       
 myCols <- cbind(condition_color,tissue_color )
 colnames(myCols)[1] <- "Condtion"
 colnames(myCols)[2] <- "Tissue"
 
 #install.packages("heatmap.plus")
 library("heatmap.plus")
 
tiff("IBN00-condition_0.01.plus.tif",width=1000,height=2000)
 heatmap.plus(hmDat, 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant",
          ColSideColors=myCols,
          scale="row",
          margins = c(10,20)) 
legend(0.8,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
legend(0.8,0.7,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F)   
dev.off()



trueR <- resCond_vV.0.01

hmDat <- logNormalizedCounts[trueR$GeneID,]
rownames(hmDat) <- trueR$Symbol

 myCols <- cbind(condition_color,tissue_color )
 colnames(myCols)[1] <- "Condtion"
 colnames(myCols)[2] <- "Tissue"
 
 #install.packages("heatmap.plus")
 library("heatmap.plus")
 FDT


############################################################################
############################################################################
############################################################################

resCond_vV.nona<- resCond_vV[!is.na(resCond_vV$Symbol),]
resCond_vV.0.001 <- resCond_vV.nona[(resCond_vV.nona$padj < 0.001),]
resCond_vV.0.001 <- resCond_vV.0.001[!is.na(resCond_vV.0.001$GeneID),]
dim(resCond_vV.0.001) 

length(unique(resCond_vV.0.001$Symbol)) # 128

length(unique(resCond_vV.0.001[resCond_vV.0.001$log2FoldChange > 0,]$Symbol))  # up 83

length(unique(resCond_vV.0.001[resCond_vV.0.001$log2FoldChange < 0,]$Symbol))  # down 45



hmDat.001 <- logNormalizedCounts[resCond_vV.0.001$GeneID,]

tiff("IBN00-condition_0.001.plus.tif",width=1000,height=2000)
 heatmap.plus(hmDat.001, 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant",
          ColSideColors=myCols,
          scale="row",
	 cexRow = 1,
	 cexCol = 2,
          margins = c(10,20)) 
legend(0.8,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
legend(0.8,0.7,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F)   
dev.off()



trueR <- resCond_vV.0.001

hmDat.001T <- logNormalizedCounts[trueR$GeneID,]
rownames(hmDat.001T) <- trueR$Symbol

 myCols <- cbind(condition_color,tissue_color )
 colnames(myCols)[1] <- "Condtion"
 colnames(myCols)[2] <- "Tissue"
 
 #install.packages("heatmap.plus")
 library("heatmap.plus")
 
tiff("IBN00-true.001.plus.tif",width=1500,height=2000)
 heatmap.plus(hmDat.001T, 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant remove tissue effect",
          ColSideColors=myCols,
          scale="row",
	 cexRow = 1,
	 cexCol = 2,
          margins = c(10,40)) 
	legend(0.88,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
	legend(0.88,0.7,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F)   
dev.off()

############################################################################
############################################################################
############################################################################

resCond_vV.nona<- resCond_vV[!is.na(resCond_vV$Symbol),]
resCond_vV.0.0001 <- resCond_vV.nona[(resCond_vV.nona$padj < 0.0001),]
resCond_vV.0.0001 <- resCond_vV.0.0001[!is.na(resCond_vV.0.0001$GeneID),]
dim(resCond_vV.0.0001) 

length(unique(resCond_vV.0.0001$Symbol)) # 37

length(unique(resCond_vV.0.0001[resCond_vV.0.0001$log2FoldChange > 0,]$Symbol))  # up 27

length(unique(resCond_vV.0.0001[resCond_vV.0.0001$log2FoldChange < 0,]$Symbol))  # down 10



hmDat.0001 <- logNormalizedCounts[resCond_vV.0.0001$GeneID,]

tiff("IBN00-condition_0.0001.plus.tif",width=1000,height=1000)
 heatmap.plus(hmDat.0001, 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant",
          ColSideColors=myCols,
          scale="row",
	 cexRow = 2,
	 cexCol = 2,
          margins = c(10,40)) 
legend(0.8,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
legend(0.8,0.7,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F)   
dev.off()



trueR <- resCond_vV.0.0001

hmDat.0001T  <- logNormalizedCounts[trueR$GeneID,]
rownames(hmDat.0001T ) <- trueR$Symbol

 myCols <- cbind(condition_color,tissue_color )
 colnames(myCols)[1] <- "Condtion"
 colnames(myCols)[2] <- "Tissue"
 
 #install.packages("heatmap.plus")
 library("heatmap.plus")
 
tiff("IBN00-true.0001.plus.tif",width=1500,height=1000)
 heatmap.plus(hmDat.0001T, 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant remove tissue effect",
          ColSideColors=myCols,
          scale="row",
	 cexRow = 2,
	 cexCol = 2,
          margins = c(10,40)) 
	legend(0.88,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
	legend(0.88,0.7,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F)   
dev.off()

