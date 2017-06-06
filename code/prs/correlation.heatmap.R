###################################################################################################################################
# This script creates heatmaps from the correlation tables
# Functions are provided to create heatmaps of the top n diseases and for a plot of the highest value per disease & cell.
###################################################################################################################################

library(grid)    
library(pheatmap)
library(RColorBrewer)

###################################################################################################################################
# functions
###################################################################################################################################
# setting rotation of colnames for pheatmap to 45 degrees
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

#function to trim names of unnecessary information
getNamesOwn <- function(list){
  sapply(strsplit(list,split = "_2015|.txt.gz|_hg19|_2016|_2014|_2013|_2012|_2011|_2010|_recessive|_additive|_mi|_EUR|_GWAS|_MLMA|_meta|_CC", perl = T),function(x){
    x[1]
  })
}

# function which creates plots of correlation data
makePlot <- function(cytTable,ct){
  plot.new()
  
  # get cell types which are used
  cytNames <- unique(cytTable[,3])
  
  # orders cyttable on different metrics: Pval, negative/positive rho scores
  cytNegCc <- cytTable[order(cytTable[,2]),]
  cytPosCc <- cytTable[rev(order(cytTable[,2])),]
  
  # which diseases are in the top 50?
  diseases <- unique(cytNegCc[,4])
  
  # create subset out of diseases from cytTable
  topTable <- cytTable[which(cytTable[,4] %in% diseases), ]
  
  # loop through diseases and get matching rho scores
  df <- sapply(diseases,function(disease){
    s <- topTable[which(topTable[,4] == disease),2]
  })
  
  rownames(df) <- cytNames
  diseases <- getNamesOwn(diseases)
  diseases <- c(paste(diseases[1:10],"neg",sep="-"),paste(diseases[11:20],"pos",sep="-"))
  colnames(df) <- paste(diseases)
  
  maxV <- max(df) + 0.1
  minV <- min(df) - 0.1
  myBreaks = seq(from = minV, to = maxV, by = 0.3)
  myColor <- colorRampPalette(c("blue","white", "red"))(16)
  
  # create heatmap of rho scores
  pdf(file = paste("/Volumes/MacOS/500fg/geneticRiskScore/plots/correlations.",ct,".diseases.pdf",sep=""), width = 11 ,height = 8.5, onefile = F)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  par(mar=c(1,1,1,1))
  pheatmap(df,show_colnames = T, show_rownames =T, fontsize_row = 5,fontsize_col = 6,
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


uniqDisHm <- function(cytTable, fontRow,ct){
  diseases <- unique(cytTable$gwas)
  cellTypes <- unique(cytTable$celltype)
  rhoDf <- cytTable
  rhoDf <- sapply(diseases,function(d){
    pheno <- cytTable[which(cytTable$gwas == d),]
    highCtScore <- sapply(cellTypes,function(ct){
      minV <- min(pheno[which(pheno$celltype == ct),5])
      maxV <- max(pheno[which(pheno$celltype == ct),5])
      if (abs(minV) < maxV){
        highest=maxV
      } else {
        highest=minV
      }
      highest
    })
    
    
  })
  
  # write.table(sort(table(cytTable[which(cytTable$pval < 0.05), 4]),decreasing = T),
  #             file = paste("/Volumes/MacOS/500fg/plots/prs/",ct,".table.tsv", sep=""), 
  #             row.names = F,
  #             quote = F,
  #             col.names = F,
  #             sep="\t")
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(rhoDf,fontsize_col = 8 ,fontsize_row = fontRow, main = ct) 
  setHook("grid.newpage", NULL, "replace")
  grid.text("Phenotype", y=-0.04, gp=gpar(fontsize=12))
  grid.text("Correlation coefficient", x= 0.8,y=0.6, gp=gpar(fontsize=7))
  grid.text(paste(ct," cell types",sep=""), x=-0.04, rot=90, gp=gpar(fontsize=12))
}
###################################################################################################################################
# data load
###################################################################################################################################
# loading data
monoTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/monocytes.correlations.table.csv",stringsAsFactors = F)
monoTable[which(is.nan(monoTable[,2])),2]<- 0


tTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/tcell.correlations.table.csv",stringsAsFactors = F)
tTable[which(is.nan(tTable[,2])),2]<- 0

cytoTable <- read.csv(file = "/Volumes/MacOS/500fg/geneticRiskScore/data/cytokine.correlations.table.csv",stringsAsFactors = F)
cytoTable[which(is.nan(cytTable[,2])),2] <- 0

# auto immune diseases in PRS scores
grepString <- "gout|disease|Disease|Major_depression|ALS|celiac_disease|Ulcerative_Colitis|celiac_disease|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma|Eczema|Type_2_Diabetes|Chronic_Kidney_Disease"

# creating plot for monocyte derived cytokines
monoRho <- monoTable[grep(grepString, monoTable[,4], ignore.case = T), ]
monoRho[which(monoRho[,1] > 0.05),2] <- 0
monoRho[,4] <- getNamesOwn(monoRho[,4])

mono.log10P <- sign(monoRho[,2])*(-log10(monoRho[,1]))
mono.log10P[is.nan(mono.log10P)] <- 0
monoRho=cbind(monoRho,mono.log10P)
makePlot(monoRho, "monocytes") #83

# creating plot for tcell derived cytokines 
tRho <- tTable[grep(grepString, tTable[,4], ignore.case = T), ]
tRho[which(tRho[,1] > 0.05),2] <- 0
tRho[,4] <- getNamesOwn(tRho[,4])

tcell.log10P <- sign(tRho[,2])*(-log10(tRho[,1]))
tcell.log10P[is.nan(tcell.log10P)] <- 0
tRho=cbind(tRho,tcell.log10P)
makePlot(tRho, "tCells")

# creating plot for tcell subset (thypea)
ibd <- "Inflammatory_Bowel_Disease|Ulcerative_Colitis|Chrohns_disease"
tIbd <- tTable[grep(ibd, tTable[,4], ignore.case = T), ]
tHypae <- tIbd[grep("hyphae", tIbd[,3]),]
makePlot(tHypae,"thypea cytokine")
makePlot(cytTable, "cytokins") #85

##########################################################################
# creates heatmap of highest value per cell per unique(disease)
##########################################################################

uniqDisHm(monoRho,5, ct = "monocytes")
# removing alzheimers since no significiant hits
tRho <- tRho[-which(tRho[,4] == "Alzheimers_disease"),]
uniqDisHm(tRho,7,ct ="tcells")




