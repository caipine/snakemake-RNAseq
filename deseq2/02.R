load ("results/preprocessing.RData")


library(ggfortify)
#library(DESeq2)
library(tidyverse)

pdf("PCA1.pdf",width=7,height=5)
# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA


autoplot(pcDat)

# We can use colour and shape to identify the Cell Type and the Status of each ample
autoplot(pcDat,
         data = sampleinfo, 
         colour="condition", 
         shape="tissue",
         size=5)
         
tiff ("IBN00.PCA.tif",width=500,height=500) 
autoplot(pcDat,
         data = sampleinfo, 
         fill="condition", 
         shape="tissue",
         size=5) +
    scale_shape_manual(values=c(21, 22,24,25)) +
    guides(fill = guide_legend(override.aes=list(shape=22)))
dev.off()

     
autoplot(pcDat,
         data = sampleinfo, 
         colour="tissue", 
         shape=FALSE,
label.size=6)

dev.off()



countVar <- apply(rlogcounts, 1, var)
# Get the row numbers for the top 500 most variable genes
highVar <- order(countVar, decreasing=TRUE)[1:3000]
# Subset logcounts matrix
hmDat <- rlogcounts[highVar,]


library(gplots)
library(RColorBrewer)

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



batch_color <- unlist( lapply(as.factor(sampleinfo$Batch), function(x){  
                  if (grepl("MCL", x)) "purple"     
                  else if (grepl("IBN", x)) "orange"
                  else if (grepl("RNA", x)) "red"
                  }))
batch_color



tissue_color <- unlist( lapply(sampleinfo$tissue, function(x){  
                  if (grepl("Spleen&LN", x)) "red"     
                  else if (grepl("PB", x)) "green"
                  else if (grepl("Apheresis", x)) "yellow"
                  else if (grepl("Spleen", x)) "blue"
                  else if (grepl("LN", x)) "purple"
                })) 
 myCols <- cbind(condition_color,tissue_color,batch_color )
 colnames(myCols)[1] <- "Condtion"
 colnames(myCols)[2] <- "Tissue"
 colnames(myCols)[3] <- "Batch"

 myCols


library(gplots)
library(RColorBrewer)
 #install.packages("heatmap.plus")
 library("heatmap.plus")
 
 pdf("whole.Heatmap_unsuper.plus.pdf",width=10,height=10)
 heatmap.plus(hmDat, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 3000 most variable genes across samples",
          ColSideColors=myCols,
          scale="row",
          margins = c(10,20)) 
legend(0.8,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen", "LN" ),fill=c("red","green","yellow","blue","purple" ),border = F)
legend(0.8,0.6,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F) 
legend(0.8,0.9,legend=c("MCL","IBN", "RNA"),fill=c("purple","orange","red"),border = F) 


dev.off()





tiff ("whole.Heatmap_unsuper.plus.tif",width=750,height=750)
 heatmap.plus(hmDat, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 3000 most variable genes across samples",
          ColSideColors=myCols,
          scale="row",
          margins = c(10,20)) 
legend(0.8,0.8,legend=c("Spleen&LN","PB","Apheresis","Spleen"),fill=c("red","green","yellow","blue"),border = F)
legend(0.8,0.6,legend=c("Resistant","Sensitive"),fill=c("purple","orange"),border = F) 
dev.off()


#save(countdata, rlogcounts, sampleinfo, file="results/preprocessing.RData")









