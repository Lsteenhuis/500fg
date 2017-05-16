#!/usr/bin/env Rscript
" 
# Date: 02-03-2017
# Author: Lars Steenhuis
# Script which checks the overlap between DE list and the drugs 
"

setwd("/Volumes/MacOS/500fg/")
drugTable <- read.csv("data/drugs.ens.csv",stringsAsFactors = F)
deResult <- list.files("deseq_results",pattern = ".csv$",full.names = T)
geneUniverse <- 64162
drugsUniverse <- nrow(drugTable)
restUniverse <- geneUniverse - drugsUniverse

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
  colProb <-lapply(deResult,itterateGen,drugGenes=drugGenes)
  return(as.numeric(colProb))  
}

"
reads DESeq result fie.
calculates the percentage overlap by taking the length of intersection between the DESeq genes and drug table subset genes.
"
itterateGen <- function(deFile,drugGenes){
  deGenes <- read.csv(deFile,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  deGenesLength <- length(deGenes.names)
  intersectedGenes <- length(intersect(deGenes.names,drugGenes))
  probability <- phyper( intersectedGenes, drugsUniverse, restUniverse, deGenesLength )
  return(as.numeric(probability))
}


runScript <- function(perc){
  print(perc)
  testProb <- apply(drugTable,2 ,itterateDrug, perc = perc)
  testProb <- as.data.frame(testProb)
  write.csv(testProb,row.names = F,file = paste("matrices/prob.de.",perc,".csv",sep=""), quote = F)
}
# Percentage that is used to get genes from the extremes
measures <- c(0.05,0.1,0.15,0.2,0.25)
measures <- c(0.1)
lapply(measures,runScript)

# get the IT number in order, used for ordering of rows.
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


colDrug <- colnames(drugTable)
colnames(prob.005.df) <- colDrug
itNa <- paste("IT",itNames,sep="")
rownames(prob.005.df) <- itNames


