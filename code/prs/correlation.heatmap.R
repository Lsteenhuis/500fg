library(grid)    
library(pheatmap)
library(RColorBrewer)

# setting rotation of colnames for pheatmap to 45 degrees
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# loading data
monoTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/monocytes.correlations.table.csv",stringsAsFactors = F)
monoTable[which(is.nan(monoTable[,2])),2]<- 0


tTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/tcell.correlations.table.csv",stringsAsFactors = F)
tTable[which(is.nan(tTable[,2])),2]<- 0

cytTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/cytokine.correlations.table.csv",stringsAsFactors = F)
cytTable[which(is.nan(cytTable[,2])),2] <- 0

#function to trim names of unnecessary information
getNamesOwn <- function(list){
  sapply(strsplit(list,split = "_2015|.txt.gz|_hg19|_2013", perl = T),function(x){
    x[1]
  })
}

# function which creates several plots of correlation data
makePlot <- function(cytTable,ct){
  # get cell types which are used
  cytNames <- unique(cytTable[,3])
  #plot.new()
  # orders cyttable on different metrics: Pval, negative/positive rho scores
  cytPval <- cytTable[order(cytTable[,1]),][1:100,]
  cyt50NegCc <- cytTable[order(cytTable[,2]),][1:10,]
  cyt50PosCc <- cytTable[rev(order(cytTable[,2])),][1:10,]
  
  par(mar=c(5,5,5,5), oma=c(0,5,0,0))
  pdf(file = paste("/Volumes/MacOS/500fg/geneticRiskScore/plots/correlations.",ct,".diseases.pdf",sep=""), width = 11 ,height = 8.5, onefile = F)
  # creates barplots of the pvalue and the rho scores or the positve and negative top scores
  barplot(rev(-log(cytPval[,1] ,base=10)), horiz = T, names.arg = getNamesOwn(cytPval[,4]), las = 1, cex.names = 0.5, main = paste("top 50 pvalues: correlation between ",ct," and prs",sep=""),
          xlab = "-log 10 of correlation P value ")
  barplot(rev(cyt50PosCc[,2]), horiz = T, names.arg = getNamesOwn(cyt50PosCc[,4]), las = 1, cex.names = 0.55, main = paste("top possitive CC: correlation between ",ct," and prs",sep=""),
          xlab = "Rho score ")
  barplot(rev(cyt50NegCc[,2]), horiz = T, names.arg = getNamesOwn(cyt50NegCc[,4]), las = 1, cex.names = 0.55, main = paste("top negative CC: correlation between ", ct," and prs",sep=""),
          xlab = "Rho score")
  #par(mar=c(5,11,5,5))
  # which diseases are in the top 50?
  diseases <- c(cyt50NegCc[,4],cyt50PosCc[,4])
  # create subset out of diseases from cytTable
  topTable <- cytTable[which(cytTable[,4] %in% diseases), ]
  # loop through diseases and get matching rho scores
  df <- sapply(diseases,function(disease){
    s <- topTable[which(topTable[,4] == disease,2),2]
  })
  
  rownames(df) <- cytNames
  #diseases <- getNamesOwn(diseases)
  diseases <- c(paste(diseases[1:10],"neg",sep="-"),paste(diseases[11:20],"pos",sep="-"))
  colnames(df) <- paste(diseases)
  
  
  myBreaks = seq(from = min(cytTable[,2]), to = max(cytTable[,2]), by = 0.5)
  myColor <- colorRampPalette(c("blue", "red"))(16)
  
  # create heatmap of rho scores
  pdf(file = paste("/Volumes/MacOS/500fg/geneticRiskScore/plots/correlations.",ct,".diseases.pdf",sep=""), width = 11 ,height = 8.5, onefile = F)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  heat <- pheatmap(df,show_colnames = T, show_rownames =T, fontsize_row = 5,fontsize_col = 6,
                   cluster_rows = T, cluster_cols = T, main =paste("top rho scores of ",ct,sep=""),
                   xlab = "gwas traits", ylab="celltypes",
                   breaks = myBreaks,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(myBreaks))
                   )
  setHook("grid.newpage", NULL, "replace")
  grid.text("GWAS traits", y=-0.04, gp=gpar(fontsize=12))
  grid.text("Rho scores", x= 0.91,y=0.95, gp=gpar(fontsize=12))
  grid.text("Cytokine cell types", x=-0.04, rot=90, gp=gpar(fontsize=12))
  dev.off()
  }


monoSubs <- monoTable[grep("disease", monoTable[,4]), ]
monoRho <- monoSubs[order(monoSubs$rho), ]
monoRho[which(monoRho[,1] > 0.05),2] <- 0
monoRho[which(sign(monoRho[,2]) == -1),2] <- -(-log(abs(monoRho[which(sign(monoRho[,2]) == -1),2] ), base = 10))
monoRho[which(sign(monoRho[,2]) == 1),2] <- -log(monoRho[which(sign(monoRho[,2]) == 1),2],base=10)

tSubs <- tTable[grep("disease", tTable[,4]), ]
tRho <- tSubs[order(tSubs$rho), ]
tRho[which(tRho[,1] > 0.05),2] <- 0
tRho[which(sign(tRho[,2]) == -1),2] <- -(-log(abs(tRho[which(sign(tRho[,2]) == -1),2] ), base = 10))
tRho[which(sign(tRho[,2]) == 1),2] <- -log(tRho[which(sign(tRho[,2]) == 1),2],base=10)


makePlot(monoRho, "monocytes") #83
makePlot(tRho, "tCells") # 86
makePlot(cytTable, "cytokins") #85


