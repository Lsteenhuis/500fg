percentages <- c("0.01", "0.02", "0.05", "0.1", "0.15", "0.2")
sapply(percentages, function(percentage){
  fp <- paste("data/gcMAP/comparison",percentage, sep = "/")
  compFiles <- list.files(fp, pattern = "drugGene", full.names = T)
  sapply(compFiles,function(cfile){
    bla <-read.csv(cfile, st)
    fileVector <- unlist(strsplit(cfile,"\\/"))
    fileN <- fileVector[length(fileVector)]
    bla$col = NA
    bla$col <- ifelse(bla[,2] > (0.05 / 114), "red","green")
    par(mfrow=c(2,3))
    #pdf(paste("plots/gCMAP/prob_pval/",percentage,"/",fileN,".pdf", sep=""))
    tit <- unlist(strsplit(fileN, "\\."))
    tit <- paste(tit[2],tit[3], sep= " ")
    plot(bla$Probability.Padj,pch = 20, col = bla$col, main = tit )
   # dev.off()
  })
})

fp = "data/gcMAP/comparison/0.01/"
compFiles <- list.files(fp, pattern = "", full.names = T)
cfile = compFiles[1]
for(cfile in compFiles){
  fileVector <- unlist(strsplit(cfile,"\\/"))
  tit <- unlist(strsplit(fileVector[6], "\\."))
  pdf(paste("plots/gCMAP/prob_pval/test/",tit[2],".",tit[3] , ".pdf", sep=""))
  par(mfrow=c(2,3))
  for(percentage in percentages){
    fileVector[4] <- percentage
    fileP <- paste(fileVector,collapse="/")
    fileN <- fileVector[length(fileVector)]
    compMatrix <- read.csv(fileP)
    #compMatrix <- compMatrix[grep("CD8\\+", compMatrix[,1], perl = T),]
    compMatrix$col = NA
    compMatrix$col <- ifelse(compMatrix[,2] > 0.05, "red","green")
    plotTitle <- paste(tit[2],tit[3], percentage, sep= " ")
    barplot(-log(compMatrix$Probability.Padj, base = 10),pch = 20, col = compMatrix$col, 
         main = plotTitle,  cex.main= 1,
         xlab = "IT", ylab="probability padj")
  }
 dev.off()
}



fp = "data/gcMAP/comparison/0.15/"
compFiles <- list.files(fp, pattern = "clio", full.names = T)
precentage = "0.15"
for(cfile in compFiles){
  fileVector <- unlist(strsplit(cfile,"\\/"))
  tit <- unlist(strsplit(fileVector[6], "\\."))
  #pdf(paste("plots/gCMAP/wilcox_pval/",tit[2],".",tit[3] , ".pdf", sep=""))
  par(mfrow=c(1,1))
  fileVector[4] <- percentage
  fileP <- paste(fileVector,collapse="/")
  fileN <- fileVector[length(fileVector)]
  compMatrix <- read.csv(fileP,stringsAsFactors = F)
  compMatrix <- compMatrix[grep("CD8\\+ Naive", compMatrix[,1], perl = T)[1:2],]
  print(compMatrix)
  compMatrix$col = NA
  compMatrix$col <- ifelse(compMatrix[,2] > 0.05, "red","green")
  plotTitle <- paste(tit[2],tit[3],"0.15", sep= " ")
  barplot(compMatrix$Probability.Padj[1:2],pch = 20, col = compMatrix$col, 
       main = plotTitle,  cex.main= 1, 
       names.arg = blanames[1:2], horiz = T )
       
  
  #dev.off()
}
  


