###################################################################################################################################
# The purpose of this script is to create grouped barplots with the help of ggplot2. 
# Plots are created of subsets from hormone + immune modulator and plateletes
###################################################################################################################################
library(ggplot2)
library(RColorBrewer)
source("/Volumes/MacOS/500fg/code/r/prs/PlottingFunctions.r")

###################################################################################################################################
# functions
###################################################################################################################################
# calculate correlations between two inhomogeneous matrices

getCorMat <- function(xMat,yMat){
  xSamples <- colnames(xMat)
  ySamples <- colnames(yMat)
  
  intersectSample <- intersect(xSamples,ySamples)
  xMat <- xMat[,which(colnames(xMat) %in% intersectSample)]
  yMat <- yMat[,which(colnames(yMat) %in% intersectSample)]
  xMat <- xMat[,order(colnames(xMat))]
  yMat <- yMat[,order(colnames(yMat))]
  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- cor(x, y,method = "spearman",use = "complete.obs")
    })
  })
  colnames(corResMat) <- rownames(xMat)
  rownames(corResMat) <- rownames(yMat)
  corResMat
  
}

# calculate correlation p values between two inhomogeneous matrices
getCorPvalMat <- function(xMat,yMat){
  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- rcorr(x, y,type = "spearman")
      res$P[2,1]
    })
  })
  colnames(corResMat) <- rownames(xMat)
  rownames(corResMat) <- rownames(yMat)
  corResMat
}
###################################################################################################################################
# loading data
###################################################################################################################################
# loading the Poly genic risk score matrix
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))

# auto immune diseases of interest

prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]

# subset of prsImmuno for figure one of the paper
figOneMatFull <- prsImmuno[grep("Rheumatoid_Arthritis|T1D|asthma",rownames(prsImmuno), ignore.case = T),]

# reading in data from the immune modulators 
immunoMod <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/20151117resultsTotalReFormat.csv", 
                        row.names = 1, header = T,sep = ",")
immunoMod <- as.data.frame(t(immunoMod))
# data needs to be log transformed
immunoMod <- log2(immunoMod)
cytokine <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv" ,
                                     header = T, row.names = 1,
                                     stringsAsFactors = F, sep = ","))
cytokine <- cytokine[,-1]



###################################################################################################################################
# preparing data
###################################################################################################################################
# create matrices of overlapping samples from figOneMat and immunoMod
allSubsets <- intersect(colnames(immunoMod), colnames(figOneMat))
figOneMatR <- figOneMat[,allSubsets]
figOneMatR <- figOneMatR[,order(allSubsets)]
immunoMod <- immunoMod[,allSubsets]
immunoMod <- immunoMod[,order(allSubsets)]
# calculate corraltions and pvalues
imCor <- getCorMat(immunoMod,figOneMatR)
imP <- getCorPvalMat(immunoMod,figOneMatR)

# create matrices of overlapping samples from figOneMat and hormones
allSubsets <- intersect(colnames(horm), colnames(figOneMat))
figOneMatR <- figOneMat[,allSubsets]
figOneMatR <- figOneMatR[,order(allSubsets)]
horm <- horm[,allSubsets]
horm <- horm[,order(allSubsets)]
# calculate correlations and pvalues
hormCor <- getCorMat(horm,figOneMatR)
hormP <- getCorPvalMat(horm,figOneMatR)
range(hormP)

# create subsets from figOneMat and plateletes data
allSubsets <- intersect(colnames(plat), colnames(figOneMat))
figOneMatR <- figOneMat[,allSubsets]
figOneMatR <- figOneMatR[,order(allSubsets)]
plat <- plat[,allSubsets]
plat <- plat[,order(allSubsets)]
# calculate correaltions and pvalues
patCor <- getCorMat(plat,figOneMatR)
patP <- getCorPvalMat(plat,figOneMatR)
range(patP)

###################################################################################################################################
# create figure one by combining correlations from immuno modulators and hormones
df_list <- list(t(imCor),t(hormCor))
imCorData <- merge_recurse(df_list)
rownames(imCorData) <- c(colnames(imCor),colnames(hormCor))
# subset imCorData and create colnames which are more readible
imCorData <- imCorData[,grep("T1D|asthma",colnames(imCorData),ignore.case = T)]
colnames(imCorData)[grep("T1D_CC",colnames(imCorData))] <- gsub("_2015_25751624_hg19.txt.gz","",colnames(imCorData)[grep("T1D_CC",colnames(imCorData))], perl=T)
colnames(imCorData)[grep("T1D_meta",colnames(imCorData))] <- gsub("_2015_25751624_hg19.txt.gz","",colnames(imCorData)[grep("T1D_meta",colnames(imCorData))], perl=T)
colnames(imCorData)[grep("Asthma_2010_860503_fixed",colnames(imCorData))] <- gsub("Asthma_2010_860503_fixed_effects_hg18_hg19.txt.gz","Asthma_FE",colnames(imCorData)[grep("Asthma_2010_860503_fixed",colnames(imCorData))], perl=T)
colnames(imCorData)[grep("Asthma_2010_860503_random",colnames(imCorData))] <- gsub("Asthma_2010_860503_random_effects_hg18_hg19.txt.gz","Asthma_RE",colnames(imCorData)[grep("Asthma_2010_860503_random",colnames(imCorData))], perl=T)

#subsetting hormone and immunomodulator of interest
imCorSub <- imCorData[grep("Adiponectin|OHP",rownames(imCorData)),]
imCorSub <- imCorSub[,order(colnames(imCorSub))]
# reordering from large to small p values
imCorSub <- imCorSub[,c("Asthma_FE_P0.01","Asthma_FE_P0.001","Asthma_FE_P1.0E-4","Asthma_FE_P1.0E-5","Asthma_FE_P5.0E-8","Asthma_RE.txt.gz_P0.01","Asthma_RE_P0.001","Asthma_RE_P1.0E-4","Asthma_RE_P1.0E-5","Asthma_RE_P5.0E-8","T1D_CC_P0.01","T1D_CC_P0.001","T1D_CC_P1.0E-4","T1D_CC_P1.0E-5" ,"T1D_CC_P5.0E-8","T1D_meta_P0.01","T1D_meta_P0.001","T1D_meta_P1.0E-4","T1D_meta_P1.0E-5","T1D_meta_P5.0E-8")  ]

pdf("/Volumes/MacOS/500fg/plots/prs/correlations.immuno.horm.pdf")
p <- bp.grouped(t(imCorSub), legend = F)
p + scale_fill_manual(values=c(positive="darkred",negative="darkblue"))  + theme(legend.position="none")
dev.off()

###################################################################################################################################
# create figure from plateletes data
patCorData <- as.data.frame(t(patCor))
patCorData <- patCorData[,grep("T1D|Rheumatoid",colnames(patCorData),ignore.case = T)]
colnames(patCorData)[grep("T1D_CC",colnames(patCorData))] <- gsub("_2015_25751624_hg19.txt.gz","",colnames(patCorData)[grep("T1D_CC",colnames(patCorData))], perl=T)
colnames(patCorData)[grep("T1D_meta",colnames(patCorData))] <- gsub("_2015_25751624_hg19.txt.gz","",colnames(patCorData)[grep("T1D_meta",colnames(patCorData))], perl=T)
colnames(patCorData)[grep("Rheumatoid_Arthritis_2010",colnames(patCorData))] <- gsub("Rheumatoid_Arthritis_2010_20453842_hg19.txt.gz","RA_2010",colnames(patCorData)[grep("Rheumatoid_Arthritis_2010",colnames(patCorData))], perl=T)
colnames(patCorData)[grep("Rheumatoid_Arthritis_2014",colnames(patCorData))] <- gsub("Rheumatoid_Arthritis_2014_24390342_hg19.txt.gz","RA_2014",colnames(patCorData)[grep("Rheumatoid_Arthritis_2014",colnames(patCorData))], perl=T)
platPlotData <- patCorData[grep("PLT_CF_AUC|PLT_CP|PLT_BTG",rownames(patCorData),ignore.case = T),]
platPlotData <- platPlotData[,order(colnames(platPlotData))]
platPlotData <- platPlotData[,c("RA_2010_P0.01","RA_2010_P0.001","RA_2010_P1.0E-4","RA_2010_P1.0E-5","RA_2010_P5.0E-8",
                                "RA_2014_P0.01","RA_2014_P0.001","RA_2014_P1.0E-4","RA_2014_P1.0E-5","RA_2014_P5.0E-8",
                                "T1D_CC_P0.01","T1D_CC_P0.001","T1D_CC_P1.0E-4","T1D_CC_P1.0E-5" ,"T1D_CC_P5.0E-8",
                                "T1D_meta_P0.01","T1D_meta_P0.001","T1D_meta_P1.0E-4","T1D_meta_P1.0E-5","T1D_meta_P5.0E-8")  ]

pdf("/Volumes/MacOS/500fg/plots/prs/cc.plat.t1d.ra.pdf")
p <- bp.grouped(t(platPlotData), legend = F)
p + scale_fill_manual(values=c(positive="darkred",negative="darkblue"))  + theme(legend.position="none")
dev.off()

###################################################################################################################################
grepImmunoString <- "Ulcerative_Colitis|Crohns_disease|Inflammatory_Bowel_Disease"
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
cytoPrs <- getCorMat(cytokine,prsImmuno)
cytoPrs <- cytoPrs[,grep("IL22|IL7|IFNy",colnames(cytoPrs),ignore.case = T)]

pheatmap(cytoPrs)
p <- bp.grouped(t(cytoPrs), legend = F)
p + scale_fill_manual(values=c(positive="darkred",negative="darkblue"))  + theme(legend.position="none")
