# This script creates a subset of brain diseases from all PRS traits
###############################################################################################
# Creating subset for  brain diseases from PRS results
# Calculating minus log 10 pvalue * sign 
###############################################################################################
nonImmune <- prs[grep("Schizophrenia|Bipolar_disorder|ALS|PGC",rownames(prs),ignore.case = T),]
myL <- getCorMat(nonImmune,cytokine)
myP <- getCorPvalMat(nonImmune,cytokine)

brainLog10Mat <- sign(myL) * (-log10(myP))
brainLog10Mat[is.na(brainLog10Mat)] <- 0
brainLog10Mat <- brainLog10Mat[which(rownames(brainLog10Mat) %in% colnames(log10PvalMat)),]

###############################################################################################
# Plotting of sign * -log10 pvalues of non immune diseases
###############################################################################################
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(brainLog10Mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(brainLog10Mat)/paletteLength, max(brainLog10Mat), length.out=floor(paletteLength/2)))

pheatmap(t(brainLog10Mat),main = "sign(correlation coefficient)  * -log10(pvalue)", fontsize = 7,
         cluster_cols = T, cluster_rows = T, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 5, cellwidth = 5,
         treeheight_row = 0,treeheight_col = 0,
         color = myColor,
         breaks = myBreaks)
###############################################################################################
#  checking amount of signicant hits of Brain and Auto immune diseases
###############################################################################################
x <-sapply(1:ncol(brainLog10Mat),function(x){
  dis <- colnames(brainLog10Mat)[x]
  sig <- length(which(abs(brainLog10Mat[,x]) > 1.3))
  insig <- length(which(abs(brainLog10Mat[,x]) < 1.3))
  list(dis,sig,insig)
})
df <- data.frame(matrix(unlist(x), nrow=ncol(brainLog10Mat), byrow=T))
colnames(df) <- c("trait","sig","insig")
length(which(abs(brainLog10Mat) > 1.3))
length(which(abs(brainLog10Mat) < 1.3))

y <-sapply(1:nrow(log10PvalMat),function(x){
  dis <- rownames(log10PvalMat)[x]
  sig <- length(which(abs(log10PvalMat[x,]) > 1.3))
  insig <- length(which(abs(log10PvalMat[x,]) < 1.3))
  list(dis,sig,insig)
})
dfy <- data.frame(matrix(unlist(y), nrow=nrow(log10PvalMat), byrow=T))
colnames(dfy) <- c("trait","sig","insig")
length(which(abs(log10PvalMat) > 1.3))
length(which(abs(log10PvalMat) < 1.3))

###############################################################################################
#  P VALUE CHECKING  fisher + t.test
# Between autoimmune diseases and brain diseases
###############################################################################################
# brain 2625 96
# other 193 1575
sigInsig <- matrix(c(83,193,
                     2542,1382),nrow=2)
sigInsig <- as.data.frame(sigInsig)
colnames(sigInsig) <- c("significant","insignificiant")
rownames(sigInsig) <- c("brain disease","auto-immune disease")

#pvalue < 2.2e-16
fisher.test(sigInsig)

id <- abs(as.vector(t(log10PvalMat[grep("E-8",rownames(log10PvalMat),perl = T),])))
bd <- abs(as.vector(brainLog10Mat[,grep("E-8",colnames(brainLog10Mat),perl = T)]))
#vialing plot
bp <-  boxplot(id,bd)
aaaaaa <- data.frame(a=c(id,bd),b=c(rep("aid",length(id)),rep("brain",length(bd))))
bla <- list(c(id,bd),c(rep("aid",length(id)),rep("brain",length(bd))))

ggplot(bla,)
p <- ggplot(MijnLijst, aes(x=MijnLijst[,1], y=MijnLijst[,2])) + 
  geom_violin() 
myP
idc <- abs(as.vector(t(prsCytoPvalMat[grep("E-8",rownames(prsCytoPvalMat),perl = T),])))
bdP <- abs(as.vector(t(myP[,grep("E-8",colnames(myP),perl = T)])))
#pvalue < 2.2e-16
boxplot(idc,bdP) 

t.test(idc,bdP)
