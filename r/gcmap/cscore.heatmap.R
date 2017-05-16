library(pheatmap)
library(RColorBrewer)
cscoreFileList <- list.files("/Volumes/MacOS/500fg/data/gcMAP/nbt_scores/D_250/", pattern="csv$", full.names = T)

scoreMatrix <- sapply(cscoreFileList, function(scoreFile){
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
