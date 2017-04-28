suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Hmisc))
registerDoParallel(3)

# read cytokine data as df and get samples
cytokine <- as.data.frame(read.csv(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv", stringsAsFactors = F))
cytoSamples <- colnames(cytokine)[3:length(colnames(cytokine))]
colnames(cytokine)[1] <- "cytokine type"

# read PRS data as df and get samples
prs <- as.data.frame(read.delim(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                sep = "\t", header = T, stringsAsFactors = F))
prsSamples <- colnames(prs)[2:ncol(prs)]

# create subset of overlapping samples
samples <- intersect(cytoSamples,prsSamples)
prs <- prs[,c(1,which(colnames(prs) %in% samples))]
cytokine <- cytokine[,c(1,which(colnames(cytokine) %in% samples))]

# divide cytokine df into monocytes and tcells subsets
monocytes <- as.data.frame(cytokine[ grep( pattern = "IL1b|IL6|TNFA", x= cytokine[,1], perl = T) , ])
tcells <- as.data.frame(cytokine[ grep( pattern = "IL17|IFNy|IL22", x = cytokine[,1], perl = T) , ])
#rm(cytokine)

# get cell names, and dataframes with only values in which Na values are set to 0 
monoNames <- monocytes[,1]
monoVal <- monocytes[,2:ncol(monocytes)]
monoVal[is.na(monoVal)] <- 0

tNames <- tcells[,1]
tVal <- tcells[,2:ncol(tcells)]
tVal[is.na(tVal)] <- 0

cytNames <- cytokine[,1]
cytVal <- cytokine[,2:ncol(cytokine)]
cytVal[is.na(cytVal)] <- 0

prsNames <- prs[,1]
prsVal <- prs[,2:ncol(prs)]

testF <- function(trait,monocytes,prs){
  # for each GWAS trait 
  as <- apply(monocytes,1,function(cellType){
    # creates a matrix with PRS and cytokine levels based on overlapping samples
    corMatrix <- cbind.data.frame(
      unlist(trait[2:431]) ,
      unlist(cellType[2:431])
    )
    # calculate correlation between PRS score and cytokine levels
    res <- rcorr(corMatrix[,1],corMatrix[,2], type = "spearman")
    #save P value, rho value, name of the cytokine, and prs name
    c(res$P[2,1],res$r[2,1],cellType[1],trait[1])
  })
  as <- t(as)
  colnames(as) <- c("pval","rho","celltype","gwas")
  as.data.frame(as,stringsAsFactors = F)
}

print("create correlation table")
cl <- makeCluster(4)
system.time({
  # for each PRS disease 
  as <- parRapply(cl,x = prs ,testF,monocytes,prs)
})
stopCluster(cl)
correlation.df <- rbind.fill(as)


cytoNames <- cytNames
#correlation.df <- createCorTable(cytokineType = cytVal)
write.table(correlation.df, "/Volumes/MacOS/500fg/geneticRiskScore/data/cytokine.correlations.table.csv", sep="," ,row.names = F, quote = F)




# create a dataframe containing correlation values for each PRS disease and celltype
createCorTable <- function(cytokineType){
  # for each PRS disease 
  as <- foreach (prsIndex = 1:nrow(prsVal)) %dopar% {
    # for each GWAS trait 
    as <- sapply(1:nrow(cytokineType),function(cytIndex){
      # creates a matrix with PRS and cytokine levels based on overlapping samples
      corMatrix <- cbind.data.frame(
        unlist(prsVal[prsIndex,]) ,
        unlist(cytokineType[cytIndex,])
      )
      # calculate correlation between PRS score and cytokine levels
      res <- rcorr(corMatrix[,1],corMatrix[,2], type = "spearman")
      #save P value, rho value, name of the cytokine, and prs name
      c(res$P[2,1],res$r[2,1],cytoNames[cytIndex],prsNames[prsIndex])
    })
    as <- t(as)
    colnames(as) <- c("pval","rho","celltype","gwas")
    as.data.frame(as,stringsAsFactors = F)
  }  
  correlation.df <- rbind.fill(as)
}