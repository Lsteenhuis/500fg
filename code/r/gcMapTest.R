library("gCMAP")
library("pheatmap")
library("DESeq")
ibrary(grid)
setwd("/Volumes/MacOS/500fg")

drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)

geneLevelExpression <- read.delim("data/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names = 1)
immuneTraits <- read.csv("data/IRT_immuneTraits_500FG.csv", row.names=1)

kegg_ensembl <- read.delim(file="/Users/umcg-lsteenhuis/Downloads/mart_export.txt",stringsAsFactors = F)
universe <- unlist(kegg_ensembl[,1])

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

perc=0.05
subsetLength <- ceiling(nrow(drugTable) * perc )
con <- c(rep("low",percentage),rep("high",percentage))

testDeSeq <- function(IT){
  percentage= ceiling(nrow(immuneTrait.subset) / 100 * 15)
  itValues <- immuneTrait.subset[,IT]
  immuneTrait.subset.order <- order(itValues,decreasing = T)
  
  lowValues <- samples[rev(immuneTrait.subset.order)[1:percentage]]
  highValues <- samples[immuneTrait.subset.order[1:percentage]]
  
  extremeValues <- c(lowValues,highValues)
  ItMatrix <- as.data.frame(geneLevelExpression[,extremeValues])
  
  return(ItMatrix)
  #return(newCountDataSet(geneLevelExpression[,extremeValues], con ))
}

generateDrugGeneSetCollection <- function(drugColIndex){
  drugCol <- drugTable[,drugColIndex]
  drugName <- colnames(drugTable)[drugColIndex]
  drug.sorted <- order(drugCol)
   
  drugGenesHi <- drugTable[,2][drug.sorted[0:subsetLength]]
  drugGenesLo <- drugTable[,2][rev(drug.sorted)[0:subsetLength]]
  
  drugGeneSetHi <- GeneSet(unique(drugGenesHi),setIdentifier = paste(drugName,"Hi",sep="_"),
                           setName =paste(drugName,"High",sep=""))
  drugGeneSetLo <- GeneSet(unique(drugGenesLo),setIdentifier = paste(drugName,"Lo",sep="_"),
                           setName =paste(drugName,"Low",sep=""))

  return(c(drugGeneSetLo,drugGeneSetHi))
}

getFishSig <- function(drugFish,drugSet) {
  drugFish <- padj(drugFish)[names(drugSet)]
  return(drugFish)
}

getFisherScores <- function(index,set) {
  drugSet <- set[index] 
  # calculate fisher score
  fisher_result_list <- fisher_score(sets = drugSet, query = cmap, universe = universe)
  # retrieve significant fisher scores
  sigFishSamples <- lapply(fisher_result_list,getFishSig, drugSet = drugSet)
  sigFishSamples <- unlist(sigFishSamples)
  log10P <- lapply(sigFishSamples, function(x){
    log10P <- -log(x, base=10)
  })
  log10P <- unlist(log10P)
  return(log10P)
}

calculateLog10 <- function(){
  drug <- 
  drugGeneIds <- drug@geneIds
  drug <- GeneSetCollection(drug)
  print(drug)
  fisher_result_list <- fisher_score(sets = as.list(drug), query = cmap, universe = universe)
  sigFishSamples <- sapply(fisher_result_list,getFishSig)
  print(length(sigFishSamples))
  # ItMatrix <- lapply(colnames(immuneTrait.subset)[1], testDeSeq)
  # lowRow <- lapply(ItMatrix, function(IT){
  #   exprM <- as.data.frame(IT)
  #   apply(exprM,1,function(geneExprRow){
  #     rowSort <- sort(geneExprRow)
  #     logRow <- log10(mean(geneExprRow[1:15])/ mean(geneExprRow[16:30]))
  #     return(logRow)
  #   })
  #})
  #print(head(ItMatrix[[1]]))
  #return(lowRow)
}



deGenesDrug <- getFisherScores(drugGeneHigh[1:2])
ItMatrix <- lapply(drugGeneHigh[1],calculateLog10)
ItMatrix <- calculateLog10(drugGeneHigh[1])

# LOGIC START
countSet <- lapply(colnames(immuneTrait.subset), testDeSeq)
cde <- generate_gCMAP_NChannelSet(countSet,
  uids=1:114,
  sample.annotation=NULL,
  platform.annotation="Entrez",
  control_perturb_col="condition",
  control="low",
  perturb="high")

cmap <- induceCMAPCollection(cde,element = "z", lower = -2 , higher = 2)
drugGeneSetCol <- sapply(3:ncol(drugTable),generateDrugGeneSetCollection)
drugGeneLow <-GeneSetCollection(unlist(drugGeneSetCol)[seq(1,2618, by=2)])
drugGeneHigh <-GeneSetCollection(unlist(drugGeneSetCol)[seq(2,2618, by=2)])

breaksList = seq(2, 40, by = 1)

p10DrugLow <-  lapply(1:length(drugGeneLow),getFisherScores,drugGeneLow)
mDrugLow <- matrix(nrow = length(p10DrugLow), ncol = 114)
for(i in 1:length(p10DrugLow)){
  mDrugLow[i,] <- p10DrugLow[[i]]
}
pheatmap(mDrugLow,labels_col = paste("IT-",1:114,sep=""),fontsize= 8, main="-log10 P value low drug expression",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrBr"))(length(breaksList)) )


p10DrugHigh <- lapply(1:length(drugGeneHigh),getFisherScores,drugGeneHigh)
mDrugHigh <- matrix(nrow = length(p10DrugHigh), ncol = 114)
for(i in 1:length(p10DrugHigh)){
  mDrugHigh[i,] <- p10DrugHigh[[i]]
}
pheatmap(mDrugHigh,labels_col = paste("IT-",1:114,sep=""),fontsize= 8, main="-log10 P value low drug expression",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrBr"))(length(breaksList)) )
# LOGIC END


draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))




