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
itNames<-lapply(deResult,function(x){
  x <- gsub(pattern="[^0-9]+", x, replacement = "")
})
regexpr("([A-Z]\\w+)",deResult,perl=T)

# DEseq results over Drug table
###############################################################################
# test <- function(deFile){
#   print(deFile)
#   deGenes <- read.csv(deFile,stringsAsFactors = F )
#   deGenes.names <- rownames(deGenes)
#   percOverlap <- apply(drugTable[,3:ncol(drugTable)], 2 , calculatePercentage, deGenes.names=deGenes.names)
#   return(as.numeric(percOverlap))
# }
# 
# calculatePercentage <- function(drugs,deGenes.names){
#   drug.sorted <- order(drugs)
#   len <- ceiling(length(drugs) * perc )
#   drugGenes <- drugTable[,2][drug.sorted[0:len]]
#   drugGenes <- c(drugGenes,drugTable[,2][rev(drug.sorted)[0:len]])
#   percentage <- (length(intersect(deGenes.names,drugGenes)) / length(deGenes.names) * 100 )
#   return(as.numeric(percentage))
# }
# 
# perc <- 0.10
# geneOverlap <- lapply(deResult[89], test)
# head(sort(unlist(moi),decreasing = T ) ,20)
# 
# for( i in seq_along(geneOverlap)){
#   barplot(sort(unlist(geneOverlap[i])),main=i,xpd = T)
# }
#############################################################################


testDrug <- function(drugCol){
  drug.sorted <- order(drugCol)
  len <- ceiling(length(drugCol) * perc )
  drugGenes <- drugTable[,2][drug.sorted[0:len]]
  drugGenes <- c(drugGenes,drugTable[,2][rev(drug.sorted)[0:len]])
  colPerc <-lapply(deResult,testGen,drugGenes=drugGenes)
  return(as.numeric(colPerc))  
}

testGen <- function(deFile,drugGenes){
  deGenes <- read.csv(deFile,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  percentage <- (length(intersect(deGenes.names,drugGenes)) / length(deGenes.names) * 100 ) 
  return(as.numeric(percentage))
}
testPerc <- apply(drugTable[,3:4],2 ,testDrug)

geneCounts <- lapply(itName.ordered,function(results){
  file <- paste("deseq_results/deseq.IT",results,".results.csv",sep="")
  deGenes <- read.csv(file,stringsAsFactors = F )
  deGenes.names <- rownames(deGenes)
  return(length(deGenes.names))
})

testPerc <- cbind(geneCount=unlist(geneCounts), testPerc)

rownames(testPerc) <- itNames
itName.ordered <- order(as.numeric(row.names(testPerc)))
testPerc <- testPerc[ order(as.numeric(row.names(testPerc))), ]
testPerc <- cbind(ImmunoTrait = paste("IT",rownames(testPerc),sep=""), testPerc)
testPerc <- as.dataframe(testPerc)


write.csv(testPerc,row.names = F,file = "matrices/drug.de.overlap.csv", quote = F)


count <- 1
apply(testPerc[,1:2],2, function(x){
  #row.names(testPerc[which(testPerc == x)])
  barplot(sort(x),main=paste("IT", count,sep=""))
  count <<- count + 1
})


