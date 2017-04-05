#!/usr/bin/env Rscript
# Date: 28-03-2017
# Author: Lars Steenhuis
# Script which uses gCMAP to test for differential expression between high and low responding genes
# Calculates the probabilty of the DE genes being part of the high / low responding genes in drugTable
# Calculates the actual enrichment score between DE genes and high / low responding genes in drugTable

suppressMessages(library("gCMAP"))
suppressMessages(library("pheatmap"))
suppressMessages(library(grid))
suppressMessages(library(RColorBrewer))

setwd("/Volumes/MacOS/500fg")
source("code/r/enrichFunctions.R")
load("data/gcMAP/nchannelSet")

######## loading and reading data
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)
geneLevelExpression <- read.delim("data/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt", row.names = 1)
immuneTraits <- read.csv("data/IRT_immuneTraits_500FG.csv", row.names=1)
it_codes <- read.csv("data/full_code_sampleInfo.unix.csv", stringsAsFactors = F)[,2]

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

# creating list of CountSet objects from the samples in immuneTrait subset
countSetList <- lapply(colnames(immuneTrait.subset), getItExprMatrix)
names(countSetList) <- paste("IT",1:114,sep="")

# data loading is quicker -> result saved in ./data/gCMAP/nchannelset
# calculates differential expression from countSetList
cdeN <- generate_gCMAP_NChannelSet(countSetList,
                                   uids=1:114,
                                   sample.annotation=NULL,
                                   platform.annotation="Entrez",
                                   control_perturb_col="condition",
                                   control="low",
                                   perturb="high")

# data loading is quicker -> results saved in ./data/gCMAP/p10Drug(Low | High)
p10DrugLow <-  lapply(1:ncol(cde),getFisherScores,drugGeneLow)
p10DrugHigh <- lapply(1:ncol(getFisherScores),drugGeneHigh)

######################################################################################## LOGIC START
# Applying threshhold on Z column of cde (throw away each z score ( < | > )  2)
cmap <- induceCMAPCollection(cde,element = "z", lower = -2 , higher = 2)

# Generating genesets for the High and Low expresssed genes based on their Rank
drugGeneSetCol <- sapply(3:ncol(drugTable),generateDrugGeneSetCollection)
drugGeneLow <-GeneSetCollection(unlist(drugGeneSetCol)[seq(1,2618, by=2)])
drugGeneHigh <-GeneSetCollection(unlist(drugGeneSetCol)[seq(2,2618, by=2)])

# calculates fisherScores -log10 (probability) between set (drugGeneHigh / drugGeneLow) and query (each column of cde)
listOfExps <- useGcmapFunction(geneDirection = "highExpr", mode = "probability")
probFiles <- list.files()

sapply()
createHeatMap()

useGcmapFunction(geneDirection = drugGeneLow, mode = "probability")


######################################################################################## LOGIC END
