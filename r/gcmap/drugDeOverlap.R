#!/usr/bin/env Rscript
" 
# Date: 02-03-2017
# Author: Lars Steenhuis
# Script which checks the overlap between DE list and the drugs 
"

setwd("/Volumes/MacOS/500fg/")
drugTable <- read.csv("data/drugs.ens.csv",stringsAsFactors = F)
deResult <- list.files("deseq_results",pattern = ".csv$",full.names = T)


"
 takes a column from the drug table
 calculates the subset length based on the percentage and retrieves that amount from high and low responders.
 loops over DESeq results file for every column to call itterateGen to calculate the overlap percentages.
 returns this percentagte.
"
itterateDrug <- function(drugCol,perc){
  drug.sorted <- order(drugCol)
  subsetLength <- ceiling(length(drugCol) * perc )
  drugGenes <- drugTable[,2][drug.sorted[0:subsetLength]]
  drugGenes <- c(drugGenes,drugTable[,2][rev(drug.sorted)[0:subsetLength]])
  colPerc <-lapply(deResult,itterateGen,drugGenes=drugGenes)
  return(as.numeric(colPerc))  
}

"
 reads DESeq result fie.
 calculates the percentage overlap by taking the length of intersection between the DESeq genes and drug table subset genes.
"
itterateGen <- function(deFile,drugGenes){
  deGenes <- read.csv(deFile,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  #print(length(intersect(deGenes.names,drugGenes)))
  #print(length(rownames(unlist(deGenes.names))))
  percentage <- length(intersect(deGenes.names,drugGenes)) / length(deGenes.names) * 100  
  return(as.numeric(percentage))
}


runScript <- function(perc){
  print(perc)
  
  # get the IT number in order, used for ordering of rows.

  
  testPerc <- apply(drugTable,2 ,itterateDrug, perc = perc)
  # reads the gene counts from the DESeq results files
  #testPerc <- as.data.frame(testPerc)
  #write.csv(testPerc,row.names = F,file = paste("matrices/DE.over.Drug.",perc,".csv",sep=""), quote = F)
}
# Percentage that is used to get genes from the extremes
measures <- c(0.05,0.1,0.15,0.2,0.25)
lapply(measures,runScript)


itNames<-lapply(deResult,function(x){
  x <- gsub(pattern="[^0-9]+", x, replacement = "")
})
itNames <- unlist(itNames)
itName.ordered <- itNames[order(as.numeric(itNames))]

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


drugTable <- read.csv("data/drugs.ens.csv",stringsAsFactors = F)
deResult <- list.files("deseq_results",pattern = ".csv$",full.names = T)
perc=0.05
test <- deResult[2]
test.R <- read.csv(test)
deGenes.names <- rownames(test.R)
drugCol = drugTable[,3]
drug.sorted <- order(drugCol)
subsetLength <- ceiling(length(drugCol) * perc )
drugGenes <- drugTable[,2][drug.sorted[0:subsetLength]]
drugGenes <- c(drugGenes,drugTable[,2][rev(drug.sorted)[0:subsetLength]])
percentage <- length(intersect(deGenes.names,drugGenes)) / length(deGenes.names) * 100 

