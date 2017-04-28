setwd("/Volumes/MacOS/500fg")
print("loading libs")
suppressMessages(library("gCMAP"))
suppressMessages(library("pheatmap"))
suppressMessages(library("DESeq"))
suppressMessages(library(grid))
suppressMessages(library(RColorBrewer))
source("code/r/enrichFunctions.R")

print("loading data objects")

load("data/gcMAP/nchannelSet")
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)
kegg_ensembl <- read.delim(file="/Users/umcg-lsteenhuis/Downloads/mart_export.txt",stringsAsFactors = F)
universe <- unlist(kegg_ensembl[,1])


# Generating genesets for the High and Low expresssed genes based on their Rank
drugGeneSetCol <- sapply(3:ncol(drugTable),generateDrugGeneSetCollection)

drugGeneLow <-GeneSetCollection(unlist(drugGeneSetCol)[seq(1,2618, by=2)])
drugGeneHigh <-GeneSetCollection(unlist(drugGeneSetCol)[seq(2,2618, by=2)])

print("starting logic")
perc = ceiling(nrow(cde) * 0.05)
sets = c(drugGeneLow, drugGeneHigh)
deGenesDireciton=c("highExpr", "lowExpr")
modes = c("wilcox", "probability")
mode= "probability"

sapply(modes, function(mode){
  print(paste("Mode: ",mode))
  sapply(deGenesDireciton,function(geneDirection){
  
    sapply(1:2,function(setIndex){
      if (setIndex == 1) {
        set = drugGeneLow
        setString = "drugGeneLow"
      } else if (setIndex == 2) {
        set = drugGeneHigh
        setString = "drugGeneHigh"
      }
      print(paste("Current object: ", geneDirection, setString, sep = " "))
      
      exp <- sapply(1:114, function(i){
        # Retrieving log_fc changes from cde data. infinite cases and inf values are set to 0
        profile <- assayDataElement(cde[,i], "log_fc")
        profile[!complete.cases(profile)] <- 0
        profile[!is.finite(profile)] <- 0
        
        # switch case for deciding between high expressed DE genes or Low expressed DE genes
        if (geneDirection == "highExpr") {
          prof.ordered <- profile[order(profile,decreasing = T)]
          names(prof.ordered) <- rownames(profile)[order(profile,decreasing = T)]
        } else {
          prof.ordered <- profile[order(profile,decreasing = F)]
          names(prof.ordered) <- rownames(profile)[order(profile,decreasing = F)]
        } 
        query <- prof.ordered[1:perc]
        # creating query from log_fc changes in prof.ordered 
        if (mode == "wilcox"){
          exp <- wilcox_score(query, set)
        } else if (mode == "probability") {
          query <- GeneSet(names(query))
          exp <- fisher_score(query = query, sets = set, universe = universe)
        }
        
        # retrieving results from exp object results
        return(exp)
      })
      fileLocation = paste("data/gcMAP/",mode,"/", sep="")
      fileName = paste(mode,geneDirection,setString, sep=".")
      print(paste("saving file: ",fileLocation, fileName, sep = ""))
      save(exp, file= paste(fileLocation, fileName, sep = ""))
      return(NULL)
    })
  })
})
print("finished logic")





