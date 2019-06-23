library(DESeq2)
library(tidyverse)

# Read the sample information into a data frame



cts6 <- read.table("counts/all.tsv006", row.names=1, header=T)
cts7 <- read.table("counts/all.tsv007", row.names=1, header=T)
seqdata.MCL <- cts6 +cts7

seqdata.IBN  <- read.table("counts/all.tsvIBN", row.names=1, header=T)

seqdata.RNA.1  <- read.table("counts/all.tsv_578_1", row.names=1, header=T)
seqdata.RNA.2  <- read.table("counts/all.tsv_578_2", row.names=1, header=T)
seqdata.RNA.3  <- read.table("counts/all.tsv_578_3", row.names=1, header=T)
#seqdata.MCL <- cts6 +cts7
seqdata <- cbind (seqdata.IBN, seqdata.MCL, seqdata.RNA.1, seqdata.RNA.2, seqdata.RNA.3)
#seqdata <- cbind (seqdata.IBN, seqdata.RNA.1, seqdata.RNA.2, seqdata.RNA.3)

#choose batch 1 and 3
###############################################
sampleinfo <- read.delim("sampleW3.txt", stringsAsFactors=F)
sampleinfo
sampleinfo<- sampleinfo[sampleinfo$tissue  %in%  c("Apheresis","PB") & sampleinfo$Batch %in% c("MCL","IBN") ,]
sampleinfo

seqdata2 <-seqdata[,colnames(seqdata) %in% sampleinfo$samples]
dim(sampleinfo)
dim(seqdata2)
dim(seqdata)

seqdata <- seqdata2
############################################

keep <- rowSums(seqdata) > 10
countdata <- seqdata[keep,]
dim(countdata)
countdata <- as.matrix(countdata)

librarySizes <- colSums(countdata)
#sampleinfo <- read.delim("sampleW.txt", stringsAsFactors=F)

sampleinfo2<- cbind(sampleinfo,data.frame(librarySizes))


tiff("wblood-bar-size-tissue.tif",width=1000,height=1000)
ggplot(data=sampleinfo2, aes(x=samples, y=librarySizes, color=tissue, fill = tissue)) +
  geom_bar(stat="identity")+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

tiff("wbllod-bar-size-batch.tif",width=1000,height=1000)
ggplot(data=sampleinfo2, aes(x=samples, y=librarySizes, color=Batch, fill = Batch)) +
  geom_bar(stat="identity")+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

tiff("w-bar-size-condition.tif",width=1000,height=1000)
ggplot(data=sampleinfo2, aes(x=samples, y=librarySizes, color=condition, fill = condition)) +
  geom_bar(stat="identity")+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()






# Get log2 counts per million
logcounts <- log2(countdata + 1)

# make a colour vector
tissueCol <- as.numeric(factor(sampleinfo$tissue)) + 1
# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=tissueCol)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")


library(ggfortify)
rlogcounts <- rlog(countdata)
save(countdata, rlogcounts, sampleinfo, file="results/preprocessing.RData")



















 ggplot(data.framelibrarySizes)), 
          col=rev(morecols(50)),
          trace="column", 
          main="Sensetive_vs_Resistant",
          ColSideColors=myCols,
          scale="row",
          margins = c(10,20)) 



tissueCol <- as.numeric(factor(sampleinfo$tissue)) + 1
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
	col=tissueCol,
        main="Barplot of library sizes")
abline(h=20e6, lty=2)



