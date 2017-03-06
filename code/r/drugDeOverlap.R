#!/usr/bin/env Rscript
#####
# Date: 02-03-2017
# Author: Lars Steenhuis
# Script which checks the overlap between DE list and the drugs 
#####

library("ggplot2")
setwd("/Volumes/MacOS/500fg/")
drugTable <- read.csv("data/drugs.ens.csv",stringsAsFactors = F)
deResult <- list.files("deseq_results",pattern = ".csv$",full.names = T)

# get the IT number in order, used for ordering of rows.
itNames<-lapply(deResult,function(x){
  x <- gsub(pattern="[^0-9]+", x, replacement = "")
})
itNames <- unlist(itNames)

# Percentage that is used to get genes from the extremes
perc <- 0.25

# takes a column from the drug table
# calculates the subset length based on the percentage and retrieves that amount from high and low responders.
# loops over DESeq results file for every column to call itterateGen to calculate the overlap percentages.
# returns this percentagte.
itterateDrug <- function(drugCol){
  drug.sorted <- order(drugCol)
  subsetLength <- ceiling(length(drugCol) * perc )
  drugGenes <- drugTable[,2][drug.sorted[0:subsetLength]]
  drugGenes <- c(drugGenes,drugTable[,2][rev(drug.sorted)[0:subsetLength]])
  colPerc <-lapply(deResult,itterateGen,drugGenes=drugGenes)
  return(as.numeric(colPerc))  
}

# reads DESeq result fie.
# calculates the percentage overlap by taking the length of intersection between the DESeq genes and drug table subset genes.
itterateGen <- function(deFile,drugGenes){
  deGenes <- read.csv(deFile,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  percentage <- (length(intersect(deGenes.names,drugGenes)) / length(drugGenes) * 100 ) 
  return(as.numeric(percentage))
}
testPerc <- apply(drugTable,2 ,itterateDrug)


itName.ordered <- itNames[order(as.numeric(itNames))]

# reads the gene counts from the DESeq results files
geneCounts <- lapply(itName.ordered,function(results){
  file <- paste("deseq_results/deseq.IT",results,".results.csv",sep="")
  deGenes <- read.csv(file,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  return(length(deGenes.names))
})

#adds column containing gene counts
testPerc <- cbind(geneCount=unlist(geneCounts), testPerc)
# adds column containing IT id's
testPerc <- cbind(ImmunoTrait = paste("IT",itName.ordered,sep=""), testPerc)
testPerc <- as.dataframe(testPerc)


write.csv(testPerc,row.names = F,file = "matrices/drug.de.overlap.25.csv", quote = F)


apply(testPerc,1, function(x){
  barplot(sort(as.numeric(x[3:1309])),main=paste(x[1],x[2],"genes",sep=" "))
})

apply(testPerc[,3:20],2, function(x){
  barplot(sort(as.numeric(x)))
  return(NULL)
})


