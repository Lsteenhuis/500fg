#!/usr/bin/env Rscript
#####
# Date: 02-22-2017
# Author: Lars Steenhuis
# Performs DESeq2 analysis on all Immune traits from IRT_immuneTraits_500FG.csv.
# Checks for differential expression of the top and bottom percentage per immune trait
#####

library("DESeq2")
setwd("/Volumes/MacOS/500fg")
geneLevelExpression <- read.delim("data/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names = 1)
immuneTraits <- read.csv("data/IRT_immuneTraits_500FG.csv", row.names=1)

# samples names from geneLevelExpression and immunutraits are made the same
newRowId <- paste(substr(rownames(immuneTraits), 1, 1), "V", substr(rownames(immuneTraits), 2,4), sep="")
rownames(immuneTraits) <- newRowId

# Retrieves the samples from which we have transcriptional profile changes 
# and creates a subset of the immunutrait levels of these samples.
# Also create an ordered subset per column.
rnaSeqSamples <- colnames(geneLevelExpression)
immuneTraits.subset <- immuneTraits[rnaSeqSamples,]
immuneTrait.subset <- immuneTraits.subset[complete.cases(immuneTraits.subset),]
samples <- rownames(immuneTrait.subset)

percentage= ceiling(length(immuneTrait.subset[,1]) / 100 * 15)

# creates a template matrix for the sample annotation used for DESeq2
perc2 <- percentage*2
condition <- c(rep("high expression",percentage),rep("low expression",percentage))
#type <- c(rep("single-read",perc2))
expression_annotation <- matrix("NA",nrow=perc2,ncol=1)
expression_annotation[,1] <- condition
#expression_annotation[,2] <- type
colnames(expression_annotation) <- c("condition")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Function which performs DESeq2 analysis per Immuno trait.
# Takes one immune Trait and order their values.
# Gets top and bottom expresion samples of this immune trait and retrieve their gene expression levels.
# 
# Sets rownames of expression_annotation to the sample names.
# Performs DESeq2, use a FDR cutoff of 0.05 and create a subset where padj < 0.05 and is ordered on the log2FoldChange.
# writes results table to results dir.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

testDeSeq <- function(IT){
  rownames(expression_annotation) <- NULL
  
  itValues <- immuneTrait.subset[,IT]
  immuneTrait.subset.order <- order(itValues,decreasing = T)
  
  highValues <- samples[immuneTrait.subset.order[1:percentage]]
  lowValues <- samples[rev(immuneTrait.subset.order)[1:percentage]]
  extremeValues <- c(highValues,lowValues)
  ItMatrix <- geneLevelExpression[,extremeValues]
  rownames(expression_annotation) <- extremeValues
  
  dds <- DESeqDataSetFromMatrix(countData = ItMatrix,
                                colData = expression_annotation,
                                design = ~ condition)
  
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.05)
  resOrdered <- res[order(res$log2FoldChange),]
  resSig <- subset(resOrdered, padj < 0.05)
  write.table(resSig,paste("deseq_fifteen/deseq.",IT,".results.csv",sep=""),sep = ",")
}

lapply(colnames(immuneTrait.subset), testDeSeq)
# for file in ./*csv; do if [[ ! $(wc -l <$file) -ge 2 ]]; then mv $file $file.nohits; fi ; done
# Check if the file contains >1 lines, if not append .nohits behind the file
