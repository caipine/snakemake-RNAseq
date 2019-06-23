
### DESeq2
library(ggfortify)
library(DESeq2)
library(tidyverse)

load ("results/preprocessing.RData")

# first lets check that our rows and columns match
all(sampleinfo$Sample == colnames(countdata))


# create the design formula
design <- as.formula(~ condition)
# create the DESeqDataSet object


ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampleinfo,
                              design = design)

ddsObj.raw$condition                             

ddsObj.raw$condition <- relevel(ddsObj.raw$condition, ref = "Sensitive")

ddsObj.raw$condition




ddsObj <- DESeq(ddsObj.raw)

                              
# Apply normalisation to DDS object
ddsObj <- estimateSizeFactors(ddsObj)
sizeFactors(ddsObj)

library(limma)
logcounts <- log2(countdata + 1)

normalizedCounts <- counts(ddsObj, normalized=TRUE) 
logNormalizedCounts <- log2(normalizedCounts + 1)


pdf("IBN.Normalization_compare.pdf",width=7,height=5)
par(mfrow=c(1,2))
plotMA(logcounts, array = 7)
abline(h=0,col="grey")
plotMA(logcounts, array = 11)
abline(h=0,col="grey")


par(mfrow=c(1,2))
plotMA(logNormalizedCounts, array = 7)
abline(h=0,col="grey")
plotMA(logNormalizedCounts, array = 11)
abline(h=0,col="grey")

dev.off()

save(countdata, rlogcounts, sampleinfo, logNormalizedCounts, ddsObj, ddsObj.raw, file="results/preprocessing_2.RData")



