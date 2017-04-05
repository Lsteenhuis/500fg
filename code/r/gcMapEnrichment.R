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
drugTable <- drugTable[isUnique(drugTable[,2]), ]
kegg_ensembl <- read.delim(file="/Users/umcg-lsteenhuis/Downloads/mart_export.txt",stringsAsFactors = F)
it_codes <- read.csv("data/full_code_sampleInfo.unix.csv", stringsAsFactors = F)[,2]
universe <- unlist(kegg_ensembl[,1])


# Generating genesets for the High and Low expresssed genes based on their Rank
drugGeneSetCol <- sapply(3:ncol(drugTable),generateDrugGeneSetCollection)

drugGeneLow <-GeneSetCollection(unlist(drugGeneSetCol)[seq(1,2618, by=2)])
drugGeneHigh <-GeneSetCollection(unlist(drugGeneSetCol)[seq(2,2618, by=2)])

print("starting logic")

sets = c(drugGeneLow, drugGeneHigh)
deGenesDireciton=c("highExpr", "lowExpr")
percentages = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)
modes = c("wilcox", "probability")
sapply(percentages, function(perc){
  print(paste("Current percentage: ", perc, sep = ""))
  sapply(modes, function(mode){
    print(paste("Current mode: ", mode, sep=""))
    subsetLength = ceiling(nrow(cde) * perc)
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
          query <- prof.ordered[1:subsetLength]
          # creating query from log_fc changes in prof.ordered 
          if (mode == "wilcox"){
            exp <- wilcox_score(query, set)
          } else if (mode == "probability") {
            query <- GeneSet(names(query))
            exp <- fisher_score(query = query, sets = set, universe = universe)
          }  else if (mode == "parametric") {
            query <- profile
            exp <- gsealm_jg_score(query = query, sets = set)
            
          }else if (mode == "hyper") {
            universe = length(universe)
            whiteballs = length(unique(drugTable[,2]))
            blackballs = universe - whiteballs
            drawnBalls = length(query)
            exp <- sapply(set, function(drug){
              inter = length(intersect(names(query), unlist(geneIds(drug) ) ) ) 
              exp <- phyper(inter,whiteballs, blackballs,drawnBalls)
            })
          }
          # retrieving results from exp object results
          return(exp)
        })
        
        fileLocation = paste("data/gcMAP/",mode,"/",perc ,"/", sep="")
        checkDir(fileLocation)
        fileName = paste(mode,geneDirection,setString, sep=".")
        print(paste("saving file: ",fileLocation, fileName, sep = ""))
        save(exp, file= paste(fileLocation, fileName, sep = ""))
        if (mode == "probability"){
          print("creating comparison files")
          createComparisonFile(perc,geneDirection,setString)
        }
        
        return(NULL)
      })
    })
  })
})
print("finished logic")
