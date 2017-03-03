#!/usr/bin/env Rscript
#####
# Date: 02-03-2017
# Author: Lars Steenhuis
# Script which checks the overlap between DE list and the drugs 
#####

setwd("/Volumes/MacOS/500fg/")
drugTable <- read.csv("data/drugs.ens.csv",stringsAsFactors = F)
deResult <- list.files("deseq_results/",pattern = ".csv$",full.names = T)


test <- function(deFile){
  print(deFile)
  deGenes <- read.csv(deFile)
  deGenes.names <- rownames(deGenes)
  muh <- apply(drugTable[,3:ncol(drugTable)], 2 , blaat, deGenes.names)
  return(as.numeric(muh))
}

blaat <- function(drugs,blaat){
  drug.sorted <- order(drugs)
  len <- ceiling(length(drugs) * perc )
  genes <- drugTable[,2][drug.sorted[0:len]]
  genes <- c(genes,drugTable[,2][rev(drug.sorted)[0:len]])
  percentage <- length(intersect(gene,genes)) / length(gene) * 100
  return(as.numeric(percentage))
}

perc <- 0.10
moi <- lapply(deResult[0:20], test)
head(sort(unlist(moi),decreasing = T ) ,20)

for(i in 1:20){
  barplot(sort(unlist(moi[i])),main=i)
}


