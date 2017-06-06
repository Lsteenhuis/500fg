# small script to generate an intermediate heatmap of the results from the NBT c score map

library(pheatmap)
library(RColorBrewer)
cscoreFileList <- list.files("/Volumes/MacOS/500fg/data/gcMAP/nbt_scores/D_250/", pattern="csv$", full.names = T)
cscoreFileList <- cscoreFileList[-which(cscoreFileList =="/Volumes/MacOS/500fg/data/gcMAP/nbt_scores/D_250//listOfDrugFiles.csv")]
scoreMatrix <- sapply(cscoreFileList, function(scoreFile){
  print(scoreFile)
  currentF <- read.csv(scoreFile,header = T,stringsAsFactors = F)
  currentF[which(currentF[,3] > 0.05),2] <- 0
  currentF[which(sign(currentF[,2]) == -1),2] <- -(-log(abs(currentF[which(sign(currentF[,2]) == -1),2] ), base = 10))
  currentF[which(sign(currentF[,2]) == 1),2] <- -log(currentF[which(sign(currentF[,2]) == 1),2],base=10)
  cscore <- currentF[,2]
  cscore
})

myBreaks = seq(from = -0.9 ,to = 0.9 , by = 0.1)
myColor <- colorRampPalette(c("blue","black", "orange"))(19)
pheatmap(scoreMatrix,
         breaks = myBreaks, color = myColor,  show_colnames = F, show_rownames = F)
