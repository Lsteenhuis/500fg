#!/usr/bin/env Rscript
#####
# Date: 02-22-2017
# Author: Lars Steenhuis
# Creates three .csv files.
# 1. Contains top 25% and bottom 25% genelevel count of the genes in the geneLevelExpression data set.
# 2. Contains the sample names of the samples which are in 1..
# 3. A DESeq2 annotation file containing information about the data.
#####

library("data.table")
setwd("~/git/500fg/")
geneLevelExpression <- read.delim("patientData/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names = 1)

percentage= length(geneLevelExpression[1,]) / 100 * 25

# This function retrieves the top 25% and bottom 25% of the gene level expression data set.
# Both will be added to the same vector which is returned to the function call.
extemeValues <- function(geneRow){
  orderedGeneRow <- order(geneRow, decreasing = T)
  newHighValues <- geneRow[orderedGeneRow][1:percentage]
  newLowValues <- geneRow[rev(orderedGeneRow)][percentage:1]
  newValues <- c(newHighValues,newLowValues)
  return(newValues)
}
# This function retrieves the sample names beloning to the top and bottom 25% of the gene level expression data set.
# They will be added to the same vector in the same way as with.
extremeSamples <- function(geneRow){
  orderedGeneRow <- order(geneRow, decreasing = T)
  samplesHigh <- colnames(geneLevelExpression)[orderedGeneRow[1:percentage]]
  samplesLow <- colnames(geneLevelExpression)[rev(orderedGeneRow)[percentage:1]]
  newValues <- c(samplesHigh,samplesLow)
  return(newValues)
}

# Calls functions to retrieve the top 25% and bottom 25% values and their samples.
# Results are transposed (so genes are rows instead of columns) and casted to a data frame
extremeValues <- apply(geneLevelExpression,1, getExtemeValues)
extremeValues <- as.data.frame(t(extremeValues))
extremeSamples <- apply(geneLevelExpression,1, getExtremeSamples)
extremeSamples <- as.data.frame(t(extremeSamples))

# These are the column names for the extremeValues and extremeSamples.
# They are also used for the annotation file used for DESeq2.
colName <- c("high1","high2","high3","high4","high5","high6","high7","high8","high9","high10","high11","high12","high13",
             "high14","high15","high16","high17","high18","high19","high20","high21","high22","high23","high24","high25",
             "low1","low2","low3","low4","low5","low6","low7","low8","low9","low10","low11","low12","low13",
             "low14","low15","low16","low17","low18","low19","low20","low21","low22","low23","low24","low25")
colnames(extremeValues) <- colName
colnames(extremeSamples) <- colName

# extremeValues and extremeValues are written as .csv files to the disk
write.csv(extremeValues, file="matrices/cellAbundance/extremes.expression.csv",row.names = T, col.names = T)
write.csv(extremeValues, file="matrices/cellAbundance/extremes.samples.csv",row.names = T)

# condition and type are columns for the annotation file which is used for DESeq2
condition <- c(rep("high expression",25),rep("low expression",25))
type <- c(rep("single-read",50))
expression_annotation <- matrix("NA",nrow=50,ncol=2)
expression_annotation[,1] <- condition
expression_annotation[,2] <- type
rownames(expression_annotation) <- colName
colnames(expression_annotation) <- c("condition","type")
write.csv(expression_annotation, file="matrices/cellAbundance/extremes.annotation.csv",row.names = T)