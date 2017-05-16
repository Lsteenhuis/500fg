library(doParallel)
library(plyr)
library(reshape2)
setwd("/Volumes/MacOS/500fg/")
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)
drugTable <- drugTable[isUnique(drugTable[,2]), ]
drugGenes <- unique(drugTable[,2])
print(Sys.time())
registerDoParallel(3)
as <- foreach (colIndex = 3:103) %dopar% {
  scores <- sapply(1:114, function(cdeIndex){
    profile <- assayDataElement(cde[,cdeIndex], "log_fc")
    profile[!complete.cases(profile)] <- 0
    profile[!is.finite(profile)] <- 0
    
    
    rep <- c(1:25)
    subsL <- profile[order(profile)][rep]
    mL <- matrix(ncol=2,nrow=length(rep))
    mL[,1] <- rownames(profile)[order(profile)][rep]
    mL[,2] <- subsL
    colnames(mL) <- c("GeneID","value")
    
    subsH <- profile[order(profile, decreasing = T)][rep]
    mH <- matrix(ncol=2,nrow=length(rep))
    mH[,1] <- rownames(profile)[rev(order(profile))][rep]
    mH[,2] <- subsH
    colnames(mH) <- c("GeneID","value")
    
    #profile <- drugTable[,colIndex]
    fDrugMat <- matrix(ncol=2, nrow= length(profile))
    #fDrugMat[,1] <- drugTable[,2]
    fDrugMat[,1] <- rownames(profile)
    fDrugMat[,2] <- profile
    colnames(fDrugMat) <- c("GeneID","value")
    
    cmap_score(mH,mL ,as.data.frame(fDrugMat))
  })
}
print(Sys.time())

scores.25 <- scores
scores.250 <- scores 
scores.100 <- scores

mat <- cbind(scores.25,scores.250,scores.100)
pheatmap(mat, labels_row = paste("IT",1:114,sep=""), cex=0.7,
         cluster_rows = T, cluster_cols = F)
