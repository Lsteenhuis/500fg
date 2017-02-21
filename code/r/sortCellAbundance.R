#####
# Date: 02-17-2017
# Author: Lars Steenhuis
# Script to sort abundance of cell types from low to high.
#####

# data set describing Immune Response Trait (IRT).
# data has normalized IRT counts per individual.
# Columns are sample ID's
# Rows contain normalized counts per Immuno trait
immuneTrait <- read.csv("/Users/larssteenhuis/500FG/patientData/IRT_immuneTraits_500FG.csv", row.names=1, header=T)
newRowId = paste(substr(rownames(immuneTrait), 1, 1), "V", substr(rownames(immuneTrait), 2,4), sep="")
rownames(immuneTrait) <- newRowId

# Table Which contains Immunotrait information of immuneTrait dataset.
full_code_sampleInfo <- read.csv("~/500FG/patientData/full_code_sampleInfo.csv")

# Data set describing genelevel expression
# Data  has gene expression level counts per individual
# Columns are sample ID's 
# Rows are ENSEMBL ID's 
geneLevelExpression <- read.delim("/Users/larssteenhuis/500FG/patientData/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names=1)
# retrieve individual ID's from which RNA seq data is available
rnaSeqIndiv <- colnames(geneLevelExpression)

# create a subset of immuneTrait containing data of individuals which have RNA-Seq data available and removes na from df.
ImmuneTrait.Subset <-immuneTrait[rnaSeqIndiv,]
ImmuneTrait.Subset <- ImmuneTrait.Subset[complete.cases(ImmuneTrait.Subset),]
immuneTrait.Subset.order = apply(ImmuneTrait.Subset, 2, function(x){return(order(x, decreasing = T))})

highValues = NULL
highNames = NULL
lowValues = NULL
lowNames = NULL
columnNames = NULL
# Loops in range 1: number of columns of ImmuneTrait.subset
# Retrieves highest and lowest items in the immunotrait columns and adds them to a data frame
# Also creates a CSV with the names of the samples that are used.
for (i in 1:ncol(ImmuneTrait.Subset)){
  newHighValues = ImmuneTrait.Subset[immuneTrait.Subset.order[1:19,i], i]
  highRowNames <- rownames(ImmuneTrait.Subset[immuneTrait.Subset.order[1:19,i], ])
  newLowValues = ImmuneTrait.Subset[immuneTrait.Subset.order[95:77,i], i]
  lowRowNames <- rownames(ImmuneTrait.Subset[immuneTrait.Subset.order[95:77,i], ])
  highValues = cbind(highValues, newHighValues, highRowNames)
  lowValues = cbind(lowValues, newLowValues, lowRowNames)
  columnNames = c(columnNames, colnames(ImmuneTrait.Subset)[i], paste(colnames(ImmuneTrait.Subset)[i], "_sample", sep=""))
}
# give columns the name of the respective Immuno trait
colnames(lowValues) <- columnNames
colnames(highValues) <- columnNames

#write files containing values to folder
#write.csv(lowValues, file="matrices/cellAbundance/it.count.low.values.csv",row.names = F)
#write.csv(highValues, file="matrices/cellAbundance/it.count.high.values.csv",row.names = F)


# This function creates two data frames for the gene level expression data 
# intended for the gene level expression data. It sorts the data per row ( per gene)
# and takes the top 25% and bottom 25% expressed samples

highValuesGen = NULL
lowValuesGen = NULL
createExpressionMatrix <- function(x){
  orderedRow <- order(x, decreasing = T)
  length <- length(x) / 100 * 25
    
  newHighValues <- x[orderedRow][1:length]
  samplesHigh = colnames(geneLevelExpression)[orderedRow[1:length]]

  newLowValues <- x[rev(orderedRow)][1:length]
  samplesLow = colnames(geneLevelExpression)[rev(orderedRow)[1:length]]

  highValuesGen <<- cbind(highValuesGen,newHighValues,samplesHigh)
  lowValuesGen <<- cbind(lowValuesGen, newLowValues,samplesLow)
}

apply(geneLevelExpression,1, createExpressionMatrix)

ensemblList <- rownames(geneLevelExpression)
columnNames = NULL
createColNameMatrix <- function(ensId){
  columnNames <<- cbind(columnNames,ensId, paste(ensId,"-samples",sep=""))
}

lapply(ensemblList, createColNameMatrix)

rownames(highValuesGen) <- NULL
rownames(lowValuesGen) <- NULL
colnames(highValuesGen) <- columnNames
colnames(lowValuesGen) <- columnNamesLow




