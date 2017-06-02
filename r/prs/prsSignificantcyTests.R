# small script to check if values are significant (create better method if it is -> see prs_aid_vs_cytokine_pvalue_checker.R)
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(log10PvalMat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(log10PvalMat)/paletteLength, max(log10PvalMat), length.out=floor(paletteLength/2)))

nonImmune <- prs[grep("Schizophrenia|Bipolar_disorder|Eczema|ALS|PGC",rownames(prs),ignore.case = T),]
myL <- getCorMat(nonImmune,cytokine)
myP <- getCorPvalMat(nonImmune,cytokine)

myC <- sign(myL) * (-log10(myP))
myC[is.na(myC)] <- 0.5
myC <- myC[which(rownames(myC) %in% colnames(log10PvalMat)),]

pheatmap(t(myC),main = "sign(correlation coefficient)  * -log10(pvalue)", fontsize = 7,
         cluster_cols = T, cluster_rows = T, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 5, cellwidth = 5,
         treeheight_row = 0,treeheight_col = 0,
         color = myColor,
         breaks = myBreaks)

# brain 2625 96
# other 193 1575
###############################################################################################
#  P VALUE CHECKING  fisher + t.test
###############################################################################################
sigInsig <- matrix(c(96,193,
                     2529,1382),nrow=2)
sigInsig <- as.data.frame(sigInsig)
colnames(sigInsig) <- c("significant","insignificiant")
rownames(sigInsig) <- c("brain disease","auto-immune disease")
fisher.test(sigInsig)

id <- abs(as.vector(t(test)))
bd <- abs(as.vector(t(myC)))
bp <-  boxplot(id,bd)
t.test(id,bd)
