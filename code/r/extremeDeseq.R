library("DESeq2")
setwd("~/git/500fg/")

highCounts <- system.file("matrices/cellAbundance/high.expression.csv")

extremesCountData <- as.matrix(read.csv("matrices/cellAbundance/extremes.expression.csv", header=T, row.names = 1), row.names = "ensembl_id")
extremesColData <- as.matrix(read.csv("matrices/cellAbundance/extremes.annotation.csv", header=T, row.names = 1), row.names = "extremes")

dds <- DESeqDataSetFromMatrix(countData = extremesCountData,
                              colData = extremesColData,
                              design = ~ condition)

dds <- DESeq(dds)
res05 <- results(dds, alpha=0.05)
summary(res05)
qresOrdered <- res05[order(res05$padj),]


library("IHW")
resIHW <- results(dds, filterFun = ihw)
summary(resIHW)
IWHresOrdered <- resIHW[order(resIHW$padj),]
