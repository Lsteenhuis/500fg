#!/usr/bin/env Rscript
#####
# Date: 02-22-2017
# Author: Lars Steenhuis
# Creates three .csv files.
# 1. Contains top 25% and bottom 25% genelevel count of the genes in the geneLevelExpression data set.
# 2. Contains the sample names of the samples which are in 1..
# 3. A DESeq2 annotation file containing information about the data.
#####

library("DESeq2")
setwd("/Volumes/MacOS/500fg")
geneLevelExpression <- read.delim("patientData/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names = 1)
immuneTraits <- read.csv("patientData/IRT_immuneTraits_500FG.csv", row.names=1)
newRowId <- paste(substr(rownames(immuneTraits), 1, 1), "V", substr(rownames(immuneTraits), 2,4), sep="")
rownames(immuneTraits) <- newRowId
percentage= length(geneLevelExpression[1,]) / 100 * 25

rnaSeqSamples <- colnames(geneLevelExpression)
immuneTraits.subset <- immuneTraits[rnaSeqSamples,]
immuneTrait.subset <- immuneTraits.subset[complete.cases(immuneTraits.subset),]
immuneTrait.subset.order = apply(immuneTrait.subset, 2, function(x){return(order(x, decreasing = T))})

rownames(immuneTrait.subset[immuneTrait.subset.order[1:19,1], ])

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
extremeIt <- NULL
for (i in 1:ncol(immuneTrait.subset)){
  highValues <- rownames(immuneTrait.subset[immuneTrait.subset.order[1:percentage,i], ])
  lowValues <- rownames(immuneTrait.subset[rev(immuneTrait.subset.order[,i])[1:25],])
  column <- c(highValues,lowValues)
  extremeIt <- cbind(extremeIt,column)
}

condition <- c(rep("high expression",25),rep("low expression",25))
type <- c(rep("single-read",50))
expression_annotation <- matrix("NA",nrow=50,ncol=2)
expression_annotation[,1] <- condition
expression_annotation[,2] <- type

colnames(expression_annotation) <- c("condition","type")

for (i in 1:ncol(extremeIt)){
  rownames(expression_annotation) <- NULL
  ItMatrix <- geneLevelExpression[,extremeIt[1:50,i]]
  samples <- colnames(ItMatrix)
  rownames(expression_annotation) <- samples
  dds <- DESeqDataSetFromMatrix(countData = ItMatrix,
                                colData = expression_annotation,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.05)
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered, padj < 0.05)
  write.table(resSig,paste("deseq_results/deseq.IT",i,".results.csv",sep=""),sep = ",")
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
