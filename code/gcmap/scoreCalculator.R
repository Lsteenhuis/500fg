#############################
# this script calculates the connectivity scores between drugs and genes
#############################
library(doParallel)
library(plyr)
library(reshape2)
library(gCMAP)
############################
# load data
############################
setwd("/Volumes/MacOS/500fg/")
load("data/gcMAP/nchannelSet")
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F,header = T)
drugGenes <- unique(drugTable[,2])
shared <- intersect(drugGenes,featureNames(cde))
drugTable <- drugTable[which(shared %in% drugTable[,2]),]
drugnames <- colnames(drugTable)
drugTable <- rbind(drugnames,drugTable)
#############################
# functions
#############################

getCScore <-function(drugCol,shared){
  library(gCMAP)
  source("data/nbt.3367-S3/CMapFxns.R")
  load("data/gcMAP/nchannelSet")
  drugName = drugCol[1]
  i = drugCol[-1]
  deg_set_len = 250
  drug_signature <- matrix(nrow=length(shared),ncol =2)
  colnames(drug_signature) <- c("GeneID","value")
  drug_signature[,1] <- shared
  drug_signature[,2] <- i
  drug_signature <- as.data.frame(drug_signature)
  
  # calculate background score distribution between gene profile of the drug
  # and a random selection of ranked genes -> (null dist for p-values)
  nTrials = 1e3
  # print (Sys.time())
  rSc = rand_cmap_scores(drug_sig=drug_signature, m_genes=shared,
                         deg_set_len=deg_set_len, nTrials=nTrials)

  cell_drug_score_list = sapply(1:114, function(ct_name) {
    # p-value for scores that are zero
    score_pvalue = NA
    profile <- assayDataElement(cde[,ct_name], "log_fc")
    profile[!complete.cases(profile)] <- 0
    profile[!is.finite(profile)] <- 0
    
    rowP <- rownames(profile)[which(shared %in% rownames(profile))]
    profile <- profile[which(shared %in% rownames(profile))]
    cell_signature <- data.frame(GeneID=rowP,value=profile,stringsAsFactors = F)
    cell_signature = cell_signature[order(cell_signature$value, 
                                          decreasing=TRUE), ]
    
    # get up and down regulated genes
    deg_set_len <- 250
    cell_up = cell_signature$GeneID[1:deg_set_len]
    gsel = (nrow(cell_signature) - deg_set_len + 1):nrow(cell_signature)
    cell_down = cell_signature$GeneID[gsel]
    
    # calculate score
    score = cmap_score(cell_up, cell_down, drug_signature)
    if (score != 0) {
      score_pvalue = length(which(abs(rSc) >= abs(score)))  / nTrials
    }
    list(score=score, pvalue=score_pvalue)
  }, simplify=F)
  
  cell_drug_scores = do.call(rbind, cell_drug_score_list)
  ofile = paste("/Volumes/MacOS/500fg/data/gcMAP/nbt_scores/D_", deg_set_len, "/",drugName , ".csv", sep="")
  write.csv(cell_drug_scores, file=ofile)
}
# create cluster
cl <- makeCluster(4)
parCapply(cl = cl, drugTable[,3:ncol(drugTable)],getCScore,shared)
# close cluster
stopCluster(cl)
  
