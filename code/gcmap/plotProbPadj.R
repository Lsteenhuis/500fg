################################
# script to bar plot the porbability padj resulting from the comparison files
################################
# loading data
percentages <- c("0.01", "0.02", "0.05", "0.1", "0.15", "0.2")
fp = "data/gcMAP/comparison/0.01/"
compFiles <- list.files(fp, pattern = "clio", full.names = T)
cfile = compFiles[1]

################################
# creates plots to compare significant / insignificant hits. 
# red = insignificant
# green = significant
# expected to see a long green line over multiple percentages for significant hits
# creates plos for the comp files
for(cfile in compFiles){
  fileVector <- unlist(strsplit(cfile,"\\/"))
  tit <- unlist(strsplit(fileVector[6], "\\."))
  #pdf(paste("plots/gCMAP/prob_pval/barplot/",tit[2],".",tit[3] , ".pdf", sep=""))
  par(mfrow=c(2,3))
  for(percentage in percentages){
    fileVector[4] <- percentage
    fileP <- paste(fileVector,collapse="/")
    fileN <- fileVector[length(fileVector)]
    compMatrix <- read.csv(fileP)
    #compMatrix <- compMatrix[grep("CD8\\+ Naive", compMatrix[,1], perl = T)[1:2],]
    compMatrix$col = NA
    compMatrix$col <- ifelse(compMatrix[,4] > 0.05, "red","green")
    #compMatrix$Probability.Padj = -log(compMatrix$Probability.Padj, base = 10)
    plotTitle <- paste(tit[2],tit[3], percentage, sep= " ")
    barplot(compMatrix$Wilcox.Rank,pch = 20, col = compMatrix$col, 
         main = plotTitle,  cex.main= 1, names.arg = compMatrix[,1], las = 1)
  }
 #dev.off()
}



fp = "data/gcMAP/comparison/0.0/"

drugPattern = c("guan","trihex", "prop","spir","clio")
cellPattern = c("CD4\\+","CD8\\+","Neut","Mono","Neut")
precentage = "0.05"
cfile = compFiles[1]

for(i in 1:5){
  compFiles <- list.files(fp, pattern = "clio", full.names = T)
  for(cfile in compFiles){
    fileVector <- unlist(strsplit(cfile,"\\/"))
    tit <- unlist(strsplit(fileVector[6], "\\."))
    #pdf(paste("plots/gCMAP/wilcox_pval/",tit[2],".",tit[3] , ".pdf", sep=""))
    par(mfrow=c(3,2))
    #fileVector[4] <- percentage
    fileP <- paste(fileVector,collapse="/")
    fileN <- fileVector[length(fileVector)]
    compMatrix <- read.csv(fileP,stringsAsFactors = F)
    #compMatrix <- compMatrix[grep("CD8\\+ Naive", compMatrix[,1], perl = T)[1:2],]
    compMatrix$col = NA
    compMatrix$col <- ifelse(compMatrix[,4] > 0.05, "red","green")
    plotTitle <- paste(tit[2],tit[3],"0.15", sep= " ")
    compMatrix$Wilcox.Padj = -log(compMatrix$Wilcox.Padj, base = 10)
    barplot(compMatrix$Probability.Padj,pch = 20, col = compMatrix$col, 
         main = plotTitle )
         
    #dev.off()
  }
}


