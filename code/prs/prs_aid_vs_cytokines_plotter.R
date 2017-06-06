suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
############################################################################################################
# FUNCTIONS 
############################################################################################################
# function to get correlation from 2 unhomogenous matrices
# xMat: matrix A
# yMat: matrix B
# output: correlation matrix
getCorMat <- function(xMat,yMat){
  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- cor(x, y,method = "spearman",use = "complete.obs")
    })
  })
  colnames(corResMat) <- rownames(xMat)
  rownames(corResMat) <- rownames(yMat)
  corResMat
}
# function to get correlation from 2 unhomogenous matrices
# xMat: matrix A
# yMat: matrix B
# output: p value matrix
getCorPvalMat <- function(xMat,yMat){
  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- rcorr(x, y,
                 type = "spearman")
      res$P[2,1]
    })
  })
  colnames(corResMat) <- rownames(xMat)
  rownames(corResMat) <- rownames(yMat)
  corResMat
}

# input: two matrices
# output: list of input with intersecting columns
# use: make matrices overlap 
makeMatricesOverlap <- function(xMat,yMat){
  xSamples <- colnames(xMat)
  ySamples <- colnames(yMat)
  
  intersectSample <- intersect(xSamples,ySamples)
  xMat <- xMat[,which(colnames(xMat) %in% intersectSample)]
  yMat <- yMat[,which(colnames(yMat) %in% intersectSample)]
  xMat <- xMat[,order(colnames(xMat))]
  yMat <- yMat[,order(colnames(yMat))]
  list(xMat,yMat)
}

# input list of two matrices you want to check the equality of
# returns all.equal(df1, df2) output
checkEqual <- function(listOfMatrices){
  df1 <- as.data.frame(listOfMatrices[1])
  df1 <- df1[order(row.names(df1)),]
  df2 <- as.data.frame(listOfMatrices[2])
  df2 <- df2[order(row.names(df2)),]
  all.equal(df1,df2)
}


############################################################################################################
# LOADING DATA
############################################################################################################
# cytokine data
# rownames: cytokine name
# colnames: sample name
# values: normalized cytokine levels
cytokine <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv" ,
                                     header = T, row.names = 1,
                                     stringsAsFactors = F, sep = ","))
cytokine <- cytokine[,-1]

# PRS data
# rows: GWAS trait
# columns: sample names
# values: Risk score
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))

listOfMatrices <- makeMatricesOverlap(cytokine,prs)
cytokine <- as.data.frame(listOfMatrices[1])
prs <- as.data.frame(listOfMatrices[2])
rm(listOfMatrices)
prs[,1] <- gsub("-","\\.",prs[,1])

# divide cytokine df into monocytes and tcells subsets
monocytes <- as.data.frame(cytokine[ grep( pattern = "IL1b|IL6|TNFA", x= rownames(cytokine), perl = T) , ])
tcells <- as.data.frame(cytokine[ grep( pattern = "IL17|IFNy|IL22", x = rownames(cytokine), perl = T) , ])

# get IBD diseases from PRS
ibd <- "Inflammatory_Bowel_Disease|Ulcerative_Colitis|Crohns"
ibdPrs <- prs[grep(ibd,rownames(prs),ignore.case = T), ]

# get Alzheimers diseases from PRS
prsAd <- prs[grep("Alzheimers",rownames(prs),ignore.case = T), ]

# get rheumatoid arthriis from PRS
prsRa <- prs[grep("Rheumatoid",rownames(prs),ignore.case = T), ]

# get *hypae cells from tcells 
hyp <- tcells[grep("hyphae", rownames(tcells),perl=T, ignore.case = T),]


############################################################################################################
# LOGIC
############################################################################################################
# Calculating correlation between *hypea cells and PRS scores from IBD diseases 
############################################################################################################
tcellsHypCor <- getCorMat(ibdPrs,hyp)
# -0.07743134  0.15722679
range(tcellsHypCor)
pheatmap(tcellsHypCor, show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cluster_cols = F, cluster_rows = T)

# compare with halla results
hallaHypIbd <- read.table("/Volumes/MacOS/500fg/halla/t.hypae.VS.ibd/similarity_table.txt",row.names = 1,comment.char = "@", header = T)

# -0.07743134  0.15722679
range(hallaHypIbd)

listOfMatrices <- makeMatricesOverlap(hallaHypIbd,tcellsHypCor)
checkEqual(listOfMatrices)

############################################################################################################
# Calculating & checking correlation between cytokine cells and Alzheimers disease 
############################################################################################################
cytoPrsAdCor <- getCorMat(cytokine,prsAd)
# -0.1250406  0.1692463
range(cytoPrsAdCor)
pheatmap(t(cytoPrsAdCor),cluster_cols = F, cluster_rows = T, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5)

#compare with halla results
hallaAdCyt <- read.table("/Volumes/MacOS/500fg/halla/prsAd.cytokine/similarity_table.txt",row.names = 1,comment.char = "@", header = T)
# -0.1237958  0.1729176
range(hallaAdCyt)
checkEqual(listOfMatrices)

listOfMatrices <- makeMatricesOverlap(cytoPrsAdCor,hallaAdCyt)
ownCytoAd <- as.data.frame(listOfMatrices[1])
hallaCytoAd <- as.data.frame(listOfMatrices[2])
ownCytoAd <- ownCytoRa[order(rownames(ownCytoRa)),]
hallaCytoAd <- hallaCytoRa[order(rownames(hallaCytoRa)),]
pheatmap(ownCytoAd, main = "Own correlation comparison", fontsize = 7,
         cluster_cols = T, cluster_rows = F, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 10, cellwidth = 10)
pheatmap(hallaCytoAd, main = "HALLA correlation results", fontsize = 7,
         cluster_cols = T, cluster_rows = F, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 10, cellwidth = 10)



############################################################################################################
# Calculating & checking correlation between Rheumatoid Arthrits and cytokines
############################################################################################################
cytoPrsRaCor <- getCorMat(cytokine,prsRa)
# -0.1462555  0.1827397
range(cytoPrsRaCor)
pheatmap(t(cytoPrsRaCor),cluster_cols = F, cluster_rows = T, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5)

hallaRaCyt <- read.table("/Volumes/MacOS/500fg/halla/prsRa.cytokine/similarity_table.txt",row.names = 1,comment.char = "@", header = T)
cytoPrsRaCor <- cytoPrsRaCor[which(rownames(cytoPrsRaCor) %in% rownames(hallaRaCyt)), ]
# -0.1237958  0.1729176
range(hallaRaCyt)
listOfMatrices <- makeMatricesOverlap(cytoPrsRaCor,hallaRaCyt)
checkEqual(listOfMatrices)
ownCytoRa <- as.data.frame(listOfMatrices[1])
hallaCytoRa <- as.data.frame(listOfMatrices[2])
ownCytoRa <- ownCytoRa[order(rownames(ownCytoRa)),]
hallaCytoRa <- hallaCytoRa[order(rownames(hallaCytoRa)),]
pheatmap(ownCytoRa, main = "Own correlation comparison", fontsize = 7,
         cluster_cols = T, cluster_rows = F, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 10, cellwidth = 10)
pheatmap(hallaCytoRa, main = "HALLA correlation results", fontsize = 7,
         cluster_cols = T, cluster_rows = F, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellheight = 10, cellwidth = 10)
############################################################################################################
#  Calculating & checking correlation between normalized 500fg gene expression data and PRS results
############################################################################################################
gene.prs <- read.table("/Volumes/MacOS/500fg/halla/genecounts.prs.tsv",sep = "\t", header= T , comment.char = "@", row.names=1)
prs.gene <- read.table("/Volumes/MacOS/500fg/halla/prs.genecounts.tsv", sep ="\t", header= T , comment.char = "@", row.names=1)
listOfMatrices <- makeMatricesOverlap(gene.prs,prs.gene)
gene.prs <- as.data.frame(listOfMatrices[1])
prs.gene <- as.data.frame(listOfMatrices[2])
prsGeneCor <- getCorMat(gene.prs,prs.gene)
prsGenePval <- getCorPvalMat(gene.prs,prs.gene)
range(prsGeneCor)
#write.table(prsGeneCor,file="/Volumes/MacOS/500fg/data/halla.normGeneExpr.prs",quote = F,sep = ",")
pheatmap(t(prsGeneCor),cluster_cols = T, cluster_rows = T, show_rownames = F,show_colnames = F)#, fontsize_row = 5, fontsize_col = 5)
############################################################################################################
# Calculating correlation between cytokines and PRS scores from IBD diseases 
############################################################################################################
simTab <- read.table("/Volumes/MacOS/500fg/halla/prs.diseases.cyto.fdr.0.5/similarity_table.txt",header=T,row.names = 1, comment.char = "@")
grepString <- "Ulcerative_Colitis|T1D|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma"
prsS <- prs[grep(grepString,rownames(prs),ignore.case = T),]
cytokineS <- cytokine[which(rownames(cytokine) %in% colnames(simTab)),]
prsCytokinesCorMatGWPrs <- getCorMat(cytokine,prs[grep("E-8",rownames(prs)),])
prsCytokinesCorMat <- getCorMat(cytokine,prs)
grepString <- "gout|celiac|Ulcerative_Colitis|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma"



prsDisease <- getCorPvalMat(cytokineS,prsS)
prsDisease <- prsCytoPvalMat[grep(grepString,rownames(prsCytoPvalMat), ignore.case = T),]
myBreaks = seq(-0.2,0.2,by=0.01)

listOfMatrices <- makeMatricesOverlap(prsDisease,simTab)
prsDisease <- as.data.frame(listOfMatrices[1])
simTab <- as.data.frame(listOfMatrices[2])
prsDisease <- prsDisease[which(rownames(prsDisease) %in% rownames(simTab)),]

# Get the same rows/columns from prsDisease in prs and cytokine
prs <- prs[which(rownames(prs) %in% rownames(prsDisease)),]
cytokine <- cytokine[which(rownames(cytokine) %in% colnames(prsDisease)),]

# get Rho score between cytokine and prs
prsCytokinesCorMat <- getCorMat(cytokine,
                                prs[grep("E-8",rownames(prs),ignore.case = T),])

# get matrix containing correlation P values
prsCytoPvalMat <- getCorPvalMat(cytokine,prs)
log10PvalMat <- sign(prsCytokinesCorMat) *(-log10(prsCytoPvalMat)) 

# removing Ulcerative_Colitis_2011_21297633 -> not an EUR study
log10PvalMat <- log10PvalMat[-grep("Ulcerative_Colitis_2011_21297633", rownames(log10PvalMat)),]

paletteLength <- 50
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(log10PvalMat)/paletteLength, max(log10PvalMat), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("blue", "white", "red"))(length(myBreaks))

log10Subs <- log10PvalMat[grep("E-8",rownames(log10PvalMat),perl = T),]
log10Subs <- log10Subs[-c(6,8,9),]
pheatmap(t(log10Subs),main = "sign(correlation coefficient)  * -log10(pvalue)", 
         fontsize = 7,
         cluster_cols = F, cluster_rows = T, 
         show_rownames = T,show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         #cellheight = 5, cellwidth = 5,
         #annotation_col= cytokineInfo_Final[,c("cytokine","stimulation","time")],
         treeheight_row = 0,treeheight_col = 0,
         color = myColor,
         breaks = myBreaks)
