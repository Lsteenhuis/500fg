library("gCMAP")
library("pheatmap")
library("DESeq")
library(grid)
library(RColorBrewer)

setwd("/Volumes/MacOS/500fg")
source("code/r/enrichFunctions.R")

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

######################################################################################## LOGIC START
print(Sys.time())
percentage= ceiling(nrow(immuneTrait.subset) / 100 * 15)
countSetList <- lapply(colnames(immuneTrait.subset), getItExprMatrix, percentage = percentage)
names(countSetList) <- paste("IT",1:114,sep="")

cdeN <- generate_gCMAP_NChannelSet(countSetList,
                                   uids=1:114,
                                   sample.annotation=NULL,
                                   platform.annotation="Entrez",
                                   control_perturb_col="condition",
                                   control="low",
                                   perturb="high")
print(Sys.time())
cmap <- induceCMAPCollection(cdeN,element = "z", lower = -2 , higher = 2)
drugGeneSetCol <- sapply(3:ncol(drugTable),generateDrugGeneSetCollection)
drugGeneLow <-GeneSetCollection(unlist(drugGeneSetCol)[seq(1,2618, by=2)])
drugGeneHigh <-GeneSetCollection(unlist(drugGeneSetCol)[seq(2,2618, by=2)])

breaksList = seq(2, 40, by = 1)

p10DrugLow <-  lapply(1:length(drugGeneLow),getFisherScores,drugGeneLow)
mDrugLow <- matrix(nrow = length(p10DrugLow), ncol = 114)
for(i in 1:length(p10DrugLow)){
  mDrugLow[i,] <- p10DrugLow[[i]]
}
pdf("log10DrugLow.pdf")
pheatmap(mDrugLow,fontsize= 6, main="-log10 P value low drug expression",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrBr"))(length(breaksList)) )
dev.off()


p10DrugHigh <- lapply(1:length(drugGeneHigh),getFisherScores,drugGeneHigh)
mDrugHigh <- matrix(nrow = length(p10DrugHigh), ncol = 114)
for(i in 1:length(p10DrugHigh)){
  mDrugHigh[i,] <- p10DrugHigh[[i]]
}
pdf("log10DrugHigh.pdf")
pheatmap(mDrugHigh,fontsize= 6, main="-log10 P value high drug expression",
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrBr"))(length(breaksList)) )
dev.off()
######################################################################################## LOGIC END